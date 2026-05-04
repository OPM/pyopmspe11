# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT
# pylint: disable=C0302, R0912, R0913, R0914, R0915, R0917, E1102

"""Utiliy function for the grid and locations in the geological models."""

import csv
import numpy as np
from shapely.geometry import Point, Polygon
from shapely.prepared import prep
from alive_progress import alive_bar
from numpy.typing import NDArray

from pyopmspe11.config.config import Config
from pyopmspe11.utils.writefile import (
    create_corner_point_grid,
    write_keywords,
    opm_files,
    write_regular_spe11c_grid,
)


def generate_files(cfg: Config) -> None:
    """Handle the deck and input files generation"""
    polygons, facies, points = getpolygons(cfg)
    if cfg.grid == "corner-point":
        xc, zc, d_x, d_y, d_z, ycent, xmx, ymy = corner(cfg, points)
        if cfg.spe11 == "spe11a":
            fipnum, fluxnum = corner_point_handling_spe11a(
                cfg, polygons, facies, xc, zc
            )
            write_keywords(cfg, fipnum, fluxnum, xmx)
        else:
            fipnum, fluxnum, porv = corner_point_handling_spe11bc(
                cfg, polygons, facies, xc, zc, ymy, ycent, d_x, d_y, d_z
            )
            write_keywords(cfg, fipnum, fluxnum, xmx, None, porv)
    else:
        if cfg.spe11 == "spe11a":
            fipnum, fluxnum, xmx, zmz = structured_handling_spe11a(
                cfg, polygons, facies
            )
            write_keywords(cfg, fipnum, fluxnum, xmx, zmz)
        else:
            fipnum, fluxnum, porv, xmx, zmz = structured_handling_spe11bc(
                cfg, polygons, facies
            )
            write_keywords(cfg, fipnum, fluxnum, xmx, zmz, porv)
    opm_files(cfg)


def prepare_structured_grid(cfg: Config) -> tuple[NDArray, NDArray, NDArray]:
    """Set the regular grid parameters"""
    dims_x, dims_y, dims_z = cfg.dims
    no_cells_x, no_cells_y, no_cells_z = cfg.nxyz
    grid = cfg.grid
    spe11 = cfg.spe11
    lower = cfg.lower
    cfg.cut = 550.0 * dims_z / 1200 if lower else 0.0
    if grid == "cartesian":
        xmx = np.linspace(0, dims_x, no_cells_x + 1)
        ymy = np.linspace(0, dims_y, no_cells_y + 1)
        zmz = (
            np.linspace(dims_z - cfg.cut, dims_z, no_cells_z + 1)
            if lower
            else np.linspace(0, dims_z, no_cells_z + 1)
        )
    if grid == "tensor":
        len_x_n = len(cfg.x_n)
        len_y_n = len(cfg.y_n)
        len_z_n = len(cfg.z_n)
        x_idx = np.concatenate(
            [np.full(num, index) for index, num in enumerate(cfg.x_n)]
        )
        x_frac = np.concatenate([np.arange(1, num + 1) / num for num in cfg.x_n])
        xmx = np.concatenate(([0.0], (x_idx + x_frac) * dims_x / len_x_n))
        y_idx = np.concatenate(
            [np.full(num, index) for index, num in enumerate(cfg.y_n)]
        )
        y_frac = np.concatenate([np.arange(1, num + 1) / num for num in cfg.y_n])
        ymy = np.concatenate(([0.0], (y_idx + y_frac) * dims_y / len_y_n))
        z_idx = np.concatenate(
            [np.full(num, index) for index, num in enumerate(cfg.z_n)]
        )
        z_frac = np.concatenate([np.arange(1, num + 1) / num for num in cfg.z_n])
        zmz = (
            (dims_z - cfg.cut) + (z_idx + z_frac) * cfg.cut / len_z_n
            if lower
            else np.concatenate(([0.0], (z_idx + z_frac) * dims_z / len_z_n))
        )
        zmz = np.asarray(zmz)
        cfg.nxyz[0] = len(xmx) - 1
        cfg.nxyz[1] = len(ymy) - 1
        cfg.nxyz[2] = len(zmz) - 1
    if cfg.widthBuffer is not None:
        if (spe11 in ["spe11b", "spe11c"]) and 1.1 * cfg.widthBuffer < xmx[1]:
            xmx = np.insert(xmx, 1, cfg.widthBuffer)
            xmx = np.insert(xmx, len(xmx) - 1, xmx[-1] - cfg.widthBuffer)
            cfg.nxyz[0] += 2
        if spe11 == "spe11c" and 1.1 * cfg.widthBuffer < ymy[1]:
            ymy = np.insert(ymy, 1, cfg.widthBuffer)
            ymy = np.insert(ymy, len(ymy) - 1, ymy[-1] - cfg.widthBuffer)
            cfg.nxyz[1] += 2
    return xmx, ymy, zmz


def vertices_centers(
    xmx: NDArray, ymy: NDArray, zmz: NDArray
) -> tuple[NDArray, NDArray, NDArray]:
    """Get the axis centers and sizes for a regular grid"""
    return (
        (xmx[1:] + xmx[:-1]) / 2.0,
        (ymy[1:] + ymy[:-1]) / 2.0,
        (zmz[1:] + zmz[:-1]) / 2.0,
    )


def vertices_sizes(
    xmx: NDArray, ymy: NDArray, zmz: NDArray
) -> tuple[NDArray, NDArray, NDArray]:
    """Get the axis centers and sizes for a regular grid"""
    return (
        xmx[1:] - xmx[:-1],
        ymy[1:] - ymy[:-1],
        zmz[1:] - zmz[:-1],
    )


def structured_handling_spe11a(
    cfg: Config, polygons: list, facies: list
) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    """Geological positions in the tensor/cartesian grid for spe11a"""
    w = 1.0 if cfg.spe11 == "spe11a" else 1200.0 / 1.2
    ztopbot = cfg.dims[2] - 0.644 * w
    zmidbot = cfg.dims[2] - 0.265 * w
    if cfg.lower:
        lowpoly = get_lower_polygon(cfg)
    xmx, ymy, zmz = prepare_structured_grid(cfg)
    xcent, ycent, zcent = vertices_centers(xmx, ymy, zmz)
    nxz = cfg.nxyz[0] * cfg.nxyz[2]
    fluxnum = np.zeros(nxz, dtype="uint8")
    fipnum = np.zeros(nxz, dtype="uint8")
    with alive_bar(nxz, bar="fish") as bar_animation:
        for k in range(cfg.nxyz[2]):
            for i in range(cfg.nxyz[0]):
                bar_animation()
                gind = i + k * cfg.nxyz[0]
                n = 0
                order = polygon_search_order(zcent[k], zmidbot, ztopbot)
                for ind in order:
                    if polygons[ind].contains(Point(xcent[i], zcent[k])):
                        n = facies[ind]
                        break
                if cfg.lower:
                    if not lowpoly.contains(Point(xcent[i], zcent[k])):
                        n = 7
                fluxnum[gind] = n
                fipnum[gind] = boxes(
                    cfg,
                    xcent[i],
                    zcent[k],
                    i,
                    n,
                )
    sensors_structured_spe11abc(cfg, fipnum, xcent, ycent, zcent)
    wells_structured_spe11abc(cfg, xcent, ycent, zcent)
    return fipnum, fluxnum, xmx, zmz


def structured_handling_spe11bc(
    cfg: Config, polygons: list[Polygon], facies: list
) -> tuple[NDArray, NDArray, list, NDArray, NDArray]:
    """Geological positions in the tensor/cartesian grid for the spe11b/c"""
    assert cfg.widthBuffer is not None
    assert cfg.pvAdded is not None
    pv_l, porv, lowpoly = 0.0, [], Polygon()
    w = 1.0 if cfg.spe11 == "spe11a" else 1200.0 / 1.2
    dims_z = cfg.dims[2]
    ztopbot = dims_z - 0.644 * w
    zmidbot = dims_z - 0.265 * w
    if cfg.lower:
        lowpoly = get_lower_polygon(cfg)
    xmx, ymy, zmz = prepare_structured_grid(cfg)
    xcent, ycent, zcent = vertices_centers(xmx, ymy, zmz)
    dx, dy, dz = vertices_sizes(xmx, ymy, zmz)
    nxyz = cfg.nxyz[0] * cfg.nxyz[1] * cfg.nxyz[2]
    fluxnum = np.zeros(nxyz, dtype="uint8")
    fipnum = np.zeros(nxyz, dtype="uint8")
    no_cells_x = cfg.nxyz[0]
    no_cells_y = cfg.nxyz[1]
    no_cells_z = cfg.nxyz[2]
    spe11 = cfg.spe11
    z_c = map_z(cfg, ycent) if spe11 == "spe11c" else 0 * ycent
    with alive_bar(no_cells_x * no_cells_z, bar="fish") as bar_animation:
        for index_z in range(no_cells_z):
            value_z = zcent[index_z]
            order = polygon_search_order(value_z, zmidbot, ztopbot)
            for index_x in range(no_cells_x):
                gind = index_x + index_z * no_cells_x * no_cells_y
                bar_animation()
                value_x = xcent[index_x]
                point = Point(value_x, value_z)
                n = 0
                for i in order:
                    if polygons[i].contains(point):
                        n = facies[i]
                        break
                fluxnum[gind] = n
                fipnum[gind] = boxes(cfg, value_x, value_z - z_c[0], index_x, n)
                pv = (
                    cfg.rock[n - 1][1] * (cfg.pvAdded + cfg.widthBuffer) if n > 0 else 0
                )
                if cfg.lower and not lowpoly.contains(point):
                    porv.append(
                        f"PORV 0 {index_x+1} {index_x+1} 2* {index_z+1} {index_z+1} /"
                    )
                elif index_x == 0 and n not in (1, 7):
                    porv.append(
                        f"PORV {pv*dy[0]*dz[index_z]} 1 1 1 1 {index_z+1} {index_z+1} /"
                    )
                    pv_l = pv
                elif index_x == no_cells_x - 1 and n not in (1, 7):
                    porv.append(
                        f"PORV {pv*dy[0]*dz[index_z]} {no_cells_x} {no_cells_x} 1 1 "
                        f"{index_z+1} {index_z+1} /"
                    )
            base_offset = index_z * no_cells_x * no_cells_y
            src = base_offset
            for index_y in range(1, no_cells_y):
                dst = base_offset + index_y * no_cells_x
                fluxnum[dst : dst + no_cells_x] = fluxnum[src : src + no_cells_x]
                for index_x in range(no_cells_x):
                    value_x = xcent[index_x]
                    fipnum[dst + index_x] = boxes(
                        cfg,
                        value_x,
                        value_z - z_c[index_y],
                        index_x,
                        fluxnum[dst + index_x],
                    )
                    if cfg.lower and not lowpoly.contains(Point(value_x, value_z)):
                        pass
                    elif index_x == 0 and fluxnum[dst + index_x] not in (1, 7):
                        porv.append(
                            f"PORV {pv_l*dy[index_y]*dz[index_z]} 1 1 {index_y+1} "
                            f"{index_y+1} {index_z+1} {index_z+1} /"
                        )
                    elif index_x == no_cells_x - 1 and fluxnum[dst + index_x] not in (
                        1,
                        7,
                    ):
                        porv.append(
                            f"PORV {pv*dy[index_y]*dz[index_z]} {no_cells_x} {no_cells_x} "
                            f"{index_y+1} {index_y+1} {index_z+1} {index_z+1} /"
                        )
                src = dst
    sensors_structured_spe11abc(cfg, fipnum, xcent, ycent, zcent)
    wells_structured_spe11abc(cfg, xcent, ycent, zcent)
    if spe11 == "spe11c":
        if not cfg.lower:
            add_pv_fipnum_front_back(
                cfg, fipnum, fluxnum, porv, dx, dz, xcent, zcent, Polygon()
            )
        else:
            add_pv_fipnum_front_back(
                cfg,
                fipnum,
                fluxnum,
                porv,
                dx,
                dz,
                xcent,
                zcent,
                lowpoly,
            )
        write_regular_spe11c_grid(cfg, xmx, ymy, zmz)
    return fipnum, fluxnum, porv, xmx, zmz


def add_pv_fipnum_front_back(
    cfg: Config,
    fipnum: NDArray,
    fluxnum: NDArray,
    porv: list,
    d_x: NDArray,
    d_z: NDArray,
    xcent: NDArray,
    zcent: NDArray,
    lowpoly: Polygon,
) -> None:
    """Buffer pore volume and bc labels also on front and back boundaries."""
    assert cfg.pvAdded is not None
    assert cfg.widthBuffer is not None
    no_cells_x = cfg.nxyz[0]
    no_cells_y = cfg.nxyz[1]
    no_cells_z = cfg.nxyz[2]
    grid = cfg.grid
    rock = cfg.rock
    pv_added = cfg.pvAdded + cfg.widthBuffer
    for k in range(no_cells_z):
        for i in range(no_cells_x - 2):
            if grid != "corner-point":
                if cfg.lower:
                    if not lowpoly.contains(Point(xcent[i], zcent[k])):
                        continue
            ind = i + 1 + k * no_cells_x * no_cells_y
            n = fluxnum[ind]
            if n not in (0, 1, 7):
                pv = rock[n - 1][1] * pv_added
                ind_xz = i + 1 + k * no_cells_x if grid == "corner-point" else k
                pv_val = pv * d_x[i + 1] * d_z[ind_xz]
                porv.append(f"PORV {pv_val} {i+2} {i+2} 1 1 {k+1} {k+1} /")
                porv.append(
                    f"PORV {pv_val} {i+2} {i+2} {no_cells_y} {no_cells_y} {k+1} {k+1} /"
                )
            set_back_front_fipnums(cfg, fipnum, fluxnum, ind)


def set_back_front_fipnums(
    cfg: Config, fipnum: NDArray, fluxnum: NDArray, ind: int
) -> None:
    """
    For the front and back boundaries in spe11c:\n
    Box A: Fipnum 13\n
    Facie 1 and Box A: Fipnum 14\n
    Box B: Fipnum 15\n
    Facie 1 and Box B: Fipnum 16\n
    Box C: Fipnum 17\n
    Facie 1 and Box C: Fipnum 18\n
    """
    offset = cfg.nxyz[0] * (cfg.nxyz[1] - 1)
    val = fipnum[ind]
    if val == 2:
        fipnum[ind] = 13
        fipnum[ind + offset] = 13
    elif val == 5:
        fipnum[ind] = 14
        fipnum[ind + offset] = 14
    elif val == 3:
        fipnum[ind] = 15
        fipnum[ind + offset] = 15
    elif val == 6:
        fipnum[ind] = 16
        fipnum[ind + offset] = 16
    elif val == 4:
        fipnum[ind] = 17
        fipnum[ind + offset] = 17
    elif val == 12:
        fipnum[ind] = 18
        fipnum[ind + offset] = 18
    elif fluxnum[ind] == 1:
        fipnum[ind] = 10
        fipnum[ind + offset] = 10
    else:
        fipnum[ind] = 11
        fipnum[ind + offset] = 11


def polygon_search_order(z: float, zmidbot: float, ztopbot: float) -> list:
    """Speed up by giving the polygon order according to the search region"""
    # fmt: off
    if z > zmidbot:
        v = [25,26,31,29,12,20,1,8,19,21,14,2,3,0,4,5,6,9,10,11,15,16,17,18,22,23,24,27,28,30,7,13]
    elif z > ztopbot:
        v = [12,20,1,8,19,21,14,2,3,25,26,29,31,0,4,5,6,9,10,11,15,16,17,18,22,23,24,27,28,30,7,13]
    else:
        v = [0,4,5,6,9,10,11,15,16,17,18,22,23,24,27,28,30,7,13,12,20,1,8,19,21,14,2,3,25,26,29,31]
    return v
    # fmt: on


def corner_point_handling_spe11a(
    cfg: Config, polygons: list[Polygon], facies: list, xc: NDArray, zc: NDArray
) -> tuple[NDArray, NDArray]:
    """Locate the geological positions in the corner-point grid for the spe11a"""
    w = 1.0 if cfg.spe11 == "spe11a" else 1200.0 / 1.2
    dims_z = cfg.dims[2]
    ztopbot = dims_z - 0.644 * w
    zmidbot = dims_z - 0.265 * w
    no_cells_x = cfg.nxyz[0]
    sensors = cfg.sensors
    wells = cfg.wellCoord
    nxz = cfg.nxyz[0] * cfg.nxyz[2]
    fluxnum = np.zeros(nxz, dtype="uint8")
    fipnum = np.zeros(nxz, dtype="uint8")
    with alive_bar(nxz, bar="fish") as bar_animation:
        for i in range(nxz):
            bar_animation()
            i_x = i % no_cells_x
            k_z = i // no_cells_x
            gind = i_x + k_z * no_cells_x
            x_val = xc[i]
            z_val = zc[i]
            n = 0
            order = polygon_search_order(z_val, zmidbot, ztopbot)
            point = Point(x_val, z_val)
            for ind in order:
                if polygons[ind].contains(point):
                    n = facies[ind]
                    break
            fluxnum[gind] = n
            fipnum[gind] = boxes(cfg, x_val, z_val, i_x, n)
    dx = (xc - sensors[0][0]) ** 2
    dz = (zc + sensors[0][2] - dims_z) ** 2
    pop1 = int(np.argmin(dx + dz))
    fipnum[pop1] = 8
    i_x = pop1 % no_cells_x
    k_z = pop1 // no_cells_x
    cfg.sensorijk[0] = [i_x, 0, k_z]
    dx = (xc - wells[0][0]) ** 2
    dz = (zc - wells[0][2]) ** 2
    idwell1 = int(np.argmin(dx + dz))
    i_x = idwell1 % no_cells_x
    k_z = idwell1 // no_cells_x
    cfg.wellijk[0] = [i_x + 1, 1, k_z + 1]
    if not cfg.lower:
        dx = (xc - sensors[1][0]) ** 2
        dz = (zc + sensors[1][2] - dims_z) ** 2
        pop2 = int(np.argmin(dx + dz))
        fipnum[pop2] = 9
        i_x = pop2 % no_cells_x
        k_z = pop2 // no_cells_x
        cfg.sensorijk[1] = [i_x, 0, k_z]
        dx = (xc - wells[1][0]) ** 2
        dz = (zc - wells[1][2]) ** 2
        idwell2 = int(np.argmin(dx + dz))
        i_x = idwell2 % no_cells_x
        k_z = idwell2 // no_cells_x
        cfg.wellijk[1] = [i_x + 1, 1, k_z + 1]
    return fipnum, fluxnum


def corner_point_handling_spe11bc(
    cfg: Config,
    polygons: list[Polygon],
    facies: list,
    xc: NDArray,
    zc: NDArray,
    ymy: NDArray,
    ycent: NDArray,
    d_x: NDArray,
    d_y: NDArray,
    d_z: NDArray,
) -> tuple[NDArray, NDArray, list]:
    """Locate the geological positions in the corner-point grid for the spe11b/c"""
    assert cfg.pvAdded is not None
    assert cfg.widthBuffer is not None
    pv_l, porv = 0.0, []
    w = 1.0 if cfg.spe11 == "spe11a" else 1200.0 / 1.2
    dims_z = cfg.dims[2]
    ztopbot = dims_z - 0.644 * w
    zmidbot = dims_z - 0.265 * w
    no_cells_x = cfg.nxyz[0]
    no_cells_y = cfg.nxyz[1]
    sensors = cfg.sensors
    wells = cfg.wellCoord
    rock = cfg.rock
    spe11 = cfg.spe11
    pv_added = cfg.pvAdded + cfg.widthBuffer
    nxz = cfg.nxyz[0] * cfg.nxyz[2]
    fluxnum = np.zeros(nxz * cfg.nxyz[1], dtype="uint8")
    fipnum = np.zeros(nxz * cfg.nxyz[1], dtype="uint8")
    xtemp = np.empty(no_cells_x, dtype=xc[0].dtype)
    ztemp = np.empty(no_cells_x, dtype=zc[0].dtype)
    z_c = map_z(cfg, ycent) if spe11 == "spe11c" else 0 * ycent
    with alive_bar(nxz, bar="fish") as bar_animation:
        for i in range(nxz):
            bar_animation()
            i_x = i % no_cells_x
            k_z = i // no_cells_x
            gind = i_x + k_z * no_cells_x * no_cells_y
            x_val = xc[i]
            z_val = zc[i]
            xtemp[i_x] = x_val
            ztemp[i_x] = z_val
            n = 0
            order = polygon_search_order(z_val, zmidbot, ztopbot)
            point = Point(x_val, z_val)
            for ind in order:
                if polygons[ind].contains(point):
                    n = facies[ind]
                    break
            fluxnum[gind] = n
            fipnum[gind] = boxes(cfg, x_val, z_val - z_c[0], i_x, n)
            pv = rock[n - 1][1] * pv_added if n > 0 else 0
            if i_x == 0 and n not in (1, 7):
                porv.append(f"PORV {pv*d_y[0]*d_z[i]} 1 1 1 1 {k_z+1} {k_z+1} /")
                pv_l = pv
            elif i_x == no_cells_x - 1 and n not in (1, 7):
                porv.append(
                    f"PORV {pv*d_y[0]*d_z[i]} {no_cells_x} {no_cells_x} 1 1 {k_z+1} {k_z+1} /"
                )
            if i_x == no_cells_x - 1:
                base_offset = k_z * no_cells_x * no_cells_y
                src = base_offset
                for j in range(1, no_cells_y):
                    dst = base_offset + j * no_cells_x
                    fluxnum[dst : dst + no_cells_x] = fluxnum[src : src + no_cells_x]
                    for i_i in range(no_cells_x):
                        fipnum[dst + i_i] = boxes(
                            cfg,
                            xtemp[i_i],
                            ztemp[i_i] - z_c[j],
                            i_i,
                            fluxnum[dst + i_i],
                        )
                        if i_i == 0 and fluxnum[dst + i_i] not in (1, 7):
                            d_zl = d_z[k_z * no_cells_x + i_i]
                            porv.append(
                                f"PORV {pv_l*d_y[j]*d_zl} 1 1 {j+1} {j+1} {k_z+1} {k_z+1} /"
                            )
                        elif i_i == no_cells_x - 1 and fluxnum[dst + i_i] not in (1, 7):
                            porv.append(
                                f"PORV {pv*d_y[j]*d_z[i]} {no_cells_x} {no_cells_x} "
                                f"{j+1} {j+1} {k_z+1} {k_z+1} /"
                            )
                    src = dst
    if spe11 == "spe11c":
        add_pv_fipnum_front_back(
            cfg, fipnum, fluxnum, porv, d_x, d_z, np.empty(0), np.empty(0), Polygon()
        )
    dx = (xc - sensors[0][0]) ** 2
    dz = (zc + sensors[0][2] - dims_z) ** 2
    pop1 = int(np.argmin(dx + dz))
    dx = (xc - wells[0][0]) ** 2
    dz = (zc - wells[0][2]) ** 2
    well1 = int(np.argmin(dx + dz))
    dx = (xc - sensors[1][0]) ** 2
    dz = (zc + sensors[1][2] - dims_z) ** 2
    pop2 = int(np.argmin(dx + dz))
    dx = (xc - wells[1][0]) ** 2
    dz = (zc - wells[1][2]) ** 2
    well2 = int(np.argmin(dx + dz))
    locate_wells_sensors_cp_spe11bc(cfg, fipnum, zc, ymy, pop1, pop2, well1, well2)
    return fipnum, fluxnum, porv


def locate_wells_sensors_cp_spe11bc(
    cfg: Config,
    fipnum: NDArray,
    zc: NDArray,
    ymy: NDArray,
    pop1: int,
    pop2: int,
    well1: int,
    well2: int,
) -> None:
    """Find wells/sources and sensors ijk positions in the corner-point spe11bc"""
    no_cells_x = cfg.nxyz[0]
    no_cells_y = cfg.nxyz[1]
    i_x = well1 % no_cells_x
    k_z = well1 // no_cells_x
    well1ijk = [i_x, 0, k_z]
    i_x = well2 % no_cells_x
    k_z = well2 // no_cells_x
    well2ijk = [i_x, 0, k_z]
    i_x = pop1 % no_cells_x
    k_z = pop1 // no_cells_x
    cfg.sensorijk[0] = [i_x, 0, k_z]
    i_x = pop2 % no_cells_x
    k_z = pop2 // no_cells_x
    cfg.sensorijk[1] = [i_x, 0, k_z]
    cfg.wellijk[0] = [well1ijk[0] + 1, 1, well1ijk[2] + 1]
    cfg.wellijk[1] = [well2ijk[0] + 1, 1, well2ijk[2] + 1]
    if cfg.spe11 == "spe11c":
        assert cfg.wellCoordF is not None
        cfg.wellijkf[0] = [cfg.wellijk[0][0], 1, cfg.wellijk[0][2]]
        cfg.wellijkf[1] = [cfg.wellijk[1][0], 1, cfg.wellijk[1][2]]
        ycent = (np.array(ymy[1:]) + np.array(ymy[:-1])) / 2.0
        wji = int(np.argmin(np.abs(cfg.wellCoord[0][1] - ycent))) + 1
        wjf = int(np.argmin(np.abs(cfg.wellCoordF[0][1] - ycent))) + 1
        cfg.sensorijk[0][1] = int(np.argmin(np.abs(cfg.sensors[0][1] - ycent)))
        cfg.sensorijk[1][1] = int(np.argmin(np.abs(cfg.sensors[1][1] - ycent)))
        cfg.wellijk[0][1] = wji
        cfg.wellijk[1][1] = wji
        cfg.wellijkf[0][1] = wjf
        cfg.wellijkf[1][1] = wjf
        cfg.wellkh = []
        nz = cfg.nxyz[2]
        z_centers = np.empty(nz, dtype=zc[0].dtype)
        z_y = map_z(cfg, ycent)
        for k in range(cfg.nxyz[2]):
            z_centers[k] = zc[well1ijk[0] + k * no_cells_x]
        for j in range(cfg.wellijk[0][1], cfg.wellijkf[0][1] + 1):
            midpoints = z_centers - z_y[j - 1] - cfg.maxelevation
            cfg.wellkh.append(
                int(np.argmin(np.abs(cfg.wellCoord[0][2] - midpoints))) + 1
            )
    fipnum[
        cfg.sensorijk[0][0]
        + cfg.sensorijk[0][1] * no_cells_x
        + cfg.sensorijk[0][2] * no_cells_x * no_cells_y
    ] = 8
    if not cfg.lower:
        fipnum[
            cfg.sensorijk[1][0]
            + cfg.sensorijk[1][1] * no_cells_x
            + cfg.sensorijk[1][2] * no_cells_x * no_cells_y
        ] = 9


def boxes(cfg: Config, x_c: float, z_c: float, idx: int, fluxnum: float) -> int:
    """Find the global indices for the different boxes for the report data"""
    z_val = cfg.dims[2] + cfg.maxelevation - z_c
    if (
        (z_val >= cfg.boxb[0][2])
        & (z_val <= cfg.boxb[1][2])
        & (x_c >= cfg.boxb[0][0])
        & (x_c <= cfg.boxb[1][0])
    ):
        return check_facie1(fluxnum, 6, 3)
    if (
        (z_val >= cfg.boxc[0][2])
        & (z_val <= cfg.boxc[1][2])
        & (x_c >= cfg.boxc[0][0])
        & (x_c <= cfg.boxc[1][0])
    ):
        return check_facie1(fluxnum, 12, 4)
    if (
        (z_val >= cfg.boxa[0][2])
        & (z_val <= cfg.boxa[1][2])
        & (x_c >= cfg.boxa[0][0])
        & (x_c <= cfg.boxa[1][0])
    ):
        return check_facie1(fluxnum, 5, 2)
    if cfg.spe11 != "spe11a" and idx in (0, cfg.nxyz[0] - 1):
        return check_facie1(fluxnum, 10, 11)
    if fluxnum == 1:
        return 7
    return 1


def check_facie1(fluxnum: float, numa: int, numb: int) -> int:
    """Handle the overlaping with facie 1"""
    if fluxnum == 1:
        return numa
    return numb


def sensors_structured_spe11abc(
    cfg: Config, fipnum: NDArray, xcent: NDArray, ycent: NDArray, zcent: NDArray
) -> None:
    """Find the i,j,k sensor indices"""
    nx, ny = len(xcent), len(ycent)
    for n, _ in enumerate(cfg.sensors):
        if cfg.lower and n == 1:
            continue
        i = int(np.argmin(np.abs(cfg.sensors[n][0] - xcent)))
        j = (
            0
            if cfg.spe11 != "spe11c"
            else int(np.argmin(np.abs(cfg.sensors[n][1] - ycent)))
        )
        k = int(np.argmin(np.abs(cfg.dims[2] - cfg.sensors[n][2] - zcent)))
        fipnum[i + j * nx + k * nx * ny] = 8 + n
        cfg.sensorijk[n][0] = i
        cfg.sensorijk[n][1] = j
        cfg.sensorijk[n][2] = k


def wells_structured_spe11abc(
    cfg: Config, xcent: NDArray, ycent: NDArray, zcent: NDArray
) -> None:
    """Function to find the wells/sources index"""
    if cfg.spe11 != "spe11c":
        for j, _ in enumerate(cfg.wellCoord):
            if cfg.lower and j == 1:
                continue
            cfg.wellijk[j][0] = int(np.argmin(np.abs(cfg.wellCoord[j][0] - xcent))) + 1
            cfg.wellijk[j][1] = int(np.argmin(np.abs(cfg.wellCoord[j][1] - ycent))) + 1
            cfg.wellijk[j][2] = int(np.argmin(np.abs(cfg.wellCoord[j][2] - zcent))) + 1
    else:
        assert cfg.wellCoordF is not None
        assert cfg.wellkh is not None
        z_y = map_z(cfg, ycent)
        for j, _ in enumerate(cfg.wellCoord):
            if cfg.lower and j == 1:
                continue
            cfg.wellijk[j][0] = int(np.argmin(np.abs(cfg.wellCoord[j][0] - xcent))) + 1
            cfg.wellijkf[j][0] = (
                int(np.argmin(np.abs(cfg.wellCoordF[j][0] - xcent))) + 1
            )
            cfg.wellijk[j][1] = int(np.argmin(np.abs(cfg.wellCoord[j][1] - ycent))) + 1
            cfg.wellijkf[j][1] = (
                int(np.argmin(np.abs(cfg.wellCoordF[j][1] - ycent))) + 1
            )
            if j == 0:
                midpoints = zcent - z_y[cfg.wellijk[j][1] - 1]
                cfg.wellijk[j][2] = (
                    int(np.argmin(np.abs(cfg.wellCoord[j][2] - midpoints))) + 1
                )
            else:
                cfg.wellijk[j][2] = (
                    int(np.argmin(np.abs(cfg.wellCoord[j][2] - zcent))) + 1
                )
        for j in range(cfg.wellijk[0][1], cfg.wellijkf[0][1] + 1):
            midpoints = zcent - z_y[j - 1] - cfg.maxelevation
            cfg.wellkh.append(
                int(np.argmin(np.abs(cfg.wellCoord[0][2] - midpoints))) + 1
            )


def map_z(cfg: Config, ycent: NDArray) -> NDArray:
    """Mapping z for spe11c as funtion of the y coordinate"""
    assert cfg.elevation is not None
    assert cfg.backElevation is not None
    elevation = cfg.elevation
    dims_y = cfg.dims[1]
    back_elevation = cfg.backElevation
    scale = 0.5 * dims_y
    z_pos = (
        -elevation
        + elevation * (1.0 - (ycent / scale - 1.0) ** 2)
        - ycent * back_elevation / dims_y
    )
    return z_pos


def getpolygons(cfg: Config) -> tuple[list[Polygon], list, list[Point]]:
    """Function to create the polygons from the benchmark geo file"""
    polygons, curves, facies = [], [], []
    lines: list[list[int]] = []
    points: list[Point] = []
    facie = 0
    h_ref = 1 if cfg.spe11 == "spe11a" else 1200.0 / 1.2
    l_ref = 1 if cfg.spe11 == "spe11a" else 8400.0 / 2.8
    with open(f"{cfg.pat}/reference_mesh/points.geo", "r", encoding="utf8") as file:
        for row in csv.reader(file, delimiter=" "):
            if row[0] == "//" and not points:
                continue
            if row[0][:5] == "Point":
                points.append(
                    [
                        l_ref * float(row[2][1:-1]),
                        cfg.dims[2] - h_ref * float(row[3][:-1]),
                    ]
                )
    with open(
        f"{cfg.pat}/reference_mesh/facies_coordinates.geo", "r", encoding="utf8"
    ) as file:
        for row in csv.reader(file, delimiter=" "):
            if row[0] in ["//", "Include"] and not lines:
                continue
            if row[0][:4] == "Line":
                lines.append([int(row[2][1:-1]), int(row[3][:-2])])
            elif row[0] == "Curve":
                facies.append(facie)
                tmp = [int(row[3][1:-1])]
                for col in row[4:-1]:
                    tmp.append(int(col[:-1]))
                tmp.append(int(row[-1][:-2]))
                curves.append(tmp)
            else:
                facie += 1
    lines_list = lines
    points_list = points
    for curve in curves:
        values = []
        l0 = lines_list[curve[0] - 1]
        values.append([points_list[l0[0] - 1][0], points_list[l0[0] - 1][1]])
        values.append([points_list[l0[1] - 1][0], points_list[l0[1] - 1][1]])
        for line in curve[1:]:
            idx = abs(line) - 1
            ln = lines_list[idx]
            pidx = ln[0] - 1 if line < 0 else ln[1] - 1
            values.append([points_list[pidx][0], points_list[pidx][1]])
        values.append(values[0])
        polygons.append(Polygon(values))
    prepared_polygons = [prep(poly) for poly in polygons]
    return prepared_polygons, facies, points


def get_lower_polygon(cfg: Config) -> Polygon:
    """Get the polygon for the lower active cells"""
    lines: list[list[list[float]]] = []
    dims_x, dims_z = cfg.dims[0], cfg.dims[2]
    with open(
        f"{cfg.pat}/reference_mesh/horizonts_18.geo", "r", encoding="utf8"
    ) as file:
        lol = file.readlines()
    newline = False
    for row in lol:
        if row[0] == "P":
            if newline:
                lines.append([])
                newline = False
            for i, column in enumerate(row):
                if column == "{":
                    points = row[i + 1 :].split(",")
                    lines[-1].append(
                        [
                            float(points[0]) * dims_x / 2.8,
                            (1.2 - float(points[1]) - float(points[2][:-3]))
                            * dims_z
                            / 1.2,
                        ]
                    )
                    break
        else:
            newline = True
    lines[-4].append(lines[-1][-1])
    lines[-4].append(lines[-1][0])
    lines[-4].append(lines[-4][0])
    return Polygon(lines[-4])


def get_lines(cfg: Config, points: list[Point]):
    """Read the points in the z-surface lines"""
    lines: list[list[list[float]]] = []
    dims_x, dims_z = cfg.dims[0], cfg.dims[2]
    if len(cfg.z_n) == 18:
        with open(
            f"{cfg.pat}/reference_mesh/horizonts_18.geo", "r", encoding="utf8"
        ) as file:
            lol = file.readlines()
        newline = False
        for row in lol:
            if row[0] == "P":
                if newline:
                    lines.append([])
                    newline = False
                for i, column in enumerate(row):
                    if column == "{":
                        parts = row[i + 1 :].split(",")
                        lines[-1].append(
                            [
                                float(parts[0]) * dims_x / 2.8,
                                (1.2 - float(parts[1]) - float(parts[2][:-3]))
                                * dims_z
                                / 1.2,
                            ]
                        )
                        break
            else:
                newline = True
    else:
        with open(
            f"{cfg.pat}/reference_mesh/horizonts_11.geo", "r", encoding="utf8"
        ) as file:
            for line in csv.reader(file, delimiter=" "):
                if line[0][:4] == "Line":
                    if not lines[-1]:
                        idx0 = int(line[2][1:-1]) - 1
                        lines[-1].append([points[idx0][0], points[idx0][1]])
                    idx1 = int(line[3][:-2]) - 1
                    lines[-1].append([points[idx1][0], points[idx1][1]])
                if len(line) > 1:
                    if line[1] == "Horizont":
                        lines.append([])
        lines = lines[::-1]
    return lines


def corner(
    cfg: Config, points: list[Point]
) -> tuple[NDArray, NDArray, NDArray, NDArray, NDArray, NDArray, NDArray, NDArray]:
    """Create a SPE11 corner-point grid"""
    lines = get_lines(cfg, points)
    dims_x, dims_y = cfg.dims[0], cfg.dims[1]
    total = sum(cfg.x_n)
    xmx = np.zeros(total + 1, dtype=float)
    p = 1
    for i, n_x in enumerate(cfg.x_n):
        for j in range(n_x):
            xmx[p] = (i + (j + 1.0) / n_x) * dims_x / len(cfg.x_n)
            p += 1
    total = sum(cfg.y_n)
    ymy = np.zeros(total + 1, dtype=float)
    p = 1
    for i, n_y in enumerate(cfg.y_n):
        for j in range(n_y):
            ymy[p] = (i + (j + 1.0) / n_y) * dims_y / len(cfg.y_n)
            p += 1
    if cfg.spe11 in ("spe11b", "spe11c"):
        assert cfg.widthBuffer is not None
        if 1.1 * cfg.widthBuffer < xmx[1]:
            xmx = np.insert(xmx, 1, cfg.widthBuffer)
            xmx = np.insert(xmx, len(xmx) - 1, xmx[-1] - cfg.widthBuffer)
    if cfg.spe11 == "spe11c":
        assert cfg.widthBuffer is not None
        if 1.1 * cfg.widthBuffer < ymy[1]:
            ymy = np.insert(ymy, 1, cfg.widthBuffer)
            ymy = np.insert(ymy, len(ymy) - 1, ymy[-1] - cfg.widthBuffer)
    cfg.nxyz[1] = len(ymy) - 1
    shf = 15 if len(cfg.z_n) == 18 else 8
    shf = shf if cfg.lower else 0
    rmline = [i for i, value in enumerate(cfg.z_n) if value == 0]
    if rmline:
        lines.pop(rmline[-1])
        cfg.z_n.pop(rmline[0])
        cfg.z_n[rmline[0]] = 1
    nx = len(xmx)
    nz = len(lines[shf:])
    total = nx * nz
    xcoord = np.empty(total, dtype=float)
    zcoord = np.empty(total, dtype=float)
    p = 0
    for xcor in xmx:
        for lcor in lines[shf:]:
            xcoord[p] = xcor
            xs = np.array([ii[0] for ii in lcor])
            idx = int(np.argmin(np.abs(xs - xcor)))
            if lcor[idx][0] < xcor:
                zcoord[p] = lcor[idx][1] + (
                    (lcor[idx + 1][1] - lcor[idx][1])
                    / (lcor[idx + 1][0] - lcor[idx][0])
                ) * (xcor - lcor[idx][0])
            else:
                zcoord[p] = lcor[idx - 1][1] + (
                    (lcor[idx][1] - lcor[idx - 1][1])
                    / (lcor[idx][0] - lcor[idx - 1][0])
                ) * (xcor - lcor[idx - 1][0])
            p += 1
    res = zcoord[-1]
    n_z = int(np.where(zcoord == res)[0][0])
    res = xcoord[xcoord > 0][0]
    n_x = int(round(len(xcoord) / np.where(xcoord == res)[0][0])) - 1
    cfg.nxyz[0], cfg.nxyz[2] = n_x, n_z
    xcoord, zcoord, cfg.nxyz[0], cfg.nxyz[2] = refinement_z(
        xcoord, zcoord, cfg.nxyz[2], cfg.z_n[shf:]
    )
    xmx = np.array(xmx)
    ycent = 0.5 * (np.array(ymy)[1:] + np.array(ymy)[:-1])
    d_y = np.array(ymy)[1:] - np.array(ymy)[:-1]
    d_x = xmx[1:] - xmx[:-1]
    create_corner_point_grid(cfg, xcoord, ymy, zcoord)
    no_cells_x = cfg.nxyz[0]
    no_cells_z = cfg.nxyz[2]
    total = no_cells_x * no_cells_z
    xc = np.empty(total, dtype=float)
    zc = np.empty(total, dtype=float)
    d_z = np.empty(total, dtype=float)
    p = 0
    for k in range(no_cells_z):
        for i in range(no_cells_x):
            n = (i * (no_cells_z + 1)) + k
            m = ((i + 1) * (no_cells_z + 1)) + k
            poly = Polygon(
                [
                    [xcoord[n], zcoord[n]],
                    [xcoord[m], zcoord[m]],
                    [xcoord[m + 1], zcoord[m + 1]],
                    [xcoord[n + 1], zcoord[n + 1]],
                    [xcoord[n], zcoord[n]],
                ]
            )
            xc[p], zc[p] = poly.centroid.coords[0]
            d_z[p] = poly.area / (xcoord[m] - xcoord[n])
            p += 1
    return xc, zc, d_x, d_y, d_z, ycent, xmx, ymy


def refinement_z(
    xci: NDArray, zci: NDArray, ncz: int, znr: list[int]
) -> tuple[NDArray, NDArray, int, int]:
    """Refine grid vertically according to znr refinement factors."""
    stride = ncz + 1
    ncols = len(xci) // stride
    total = ncols * (1 + np.sum(znr))
    xcr = np.empty(total, dtype=xci.dtype)
    zcr = np.empty(total, dtype=zci.dtype)
    p = 0
    for j in range(ncols):
        b = j * stride
        xcr[p] = xci[b]
        zcr[p] = zci[b]
        p += 1
        for i in range(ncz):
            xi = xci[b + i]
            zi = zci[b + i]
            dx = xci[b + i + 1] - xi
            dz = zci[b + i + 1] - zi
            w = np.arange(1.0 / znr[i], 1.0 + 1.0 / znr[i], 1.0 / znr[i])
            xcr[p : p + znr[i]] = xi + dx * w
            zcr[p : p + znr[i]] = zi + dz * w
            p += znr[i]
    ncx_new = ncols - 1
    ncz_new = np.where(zcr == zcr[-1])[0][0]
    return xcr, zcr, int(ncx_new), int(ncz_new)
