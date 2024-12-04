# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT
# pylint: disable=C0302, R0914, R0915

"""
Utiliy function for the grid and locations in the geological models.
"""

import csv
import numpy as np
import pandas as pd
from shapely.geometry import Polygon
from resdata.resfile import ResdataFile
from resdata.grid import Grid

try:
    from opm.io.ecl import EGrid as OpmGrid
    from opm.io.ecl import EclFile as OpmFile
except ImportError:
    pass


def grid(dic):
    """
    Handle the different grid types (Cartesian, tensor, and corner-point grids)

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    if dic["grid"] == "corner-point":
        corner(dic)
    elif dic["grid"] == "cartesian":
        dic["dsize"] = [1.0 * dic["dims"][i] / dic["noCells"][i] for i in range(3)]
        for i, name in enumerate(["xmx", "ymy", "zmz"]):
            dic[f"{name}"] = np.linspace(0, dic["dims"][i], dic["noCells"][i] + 1)
    else:
        dic["xcor"], dic["zcor"] = [], []
        for i, (name, arr) in enumerate(
            zip(["xmx", "ymy", "zmz"], ["x_n", "y_n", "z_n"])
        ):
            dic[f"{name}"] = [0.0]
            for j, num in enumerate(dic[f"{arr}"]):
                for k in range(num):
                    dic[f"{name}"].append(
                        (j + (k + 1.0) / num) * dic["dims"][i] / len(dic[f"{arr}"])
                    )
            dic[f"{name}"] = np.array(dic[f"{name}"])
            dic["noCells"][i] = len(dic[f"{name}"]) - 1
    if dic["grid"] != "corner-point":
        if (dic["spe11"] == "spe11b" or dic["spe11"] == "spe11c") and 1.1 * dic[
            "widthBuffer"
        ] < dic["xmx"][1]:
            dic["xmx"] = np.insert(dic["xmx"], 1, dic["widthBuffer"])
            dic["xmx"] = np.insert(
                dic["xmx"], len(dic["xmx"]) - 1, dic["xmx"][-1] - dic["widthBuffer"]
            )
            dic["noCells"][0] += 2
        if dic["spe11"] == "spe11c" and 1.1 * dic["widthBuffer"] < dic["ymy"][1]:
            dic["ymy"] = np.insert(dic["ymy"], 1, dic["widthBuffer"])
            dic["ymy"] = np.insert(
                dic["ymy"], len(dic["ymy"]) - 1, dic["ymy"][-1] - dic["widthBuffer"]
            )
            dic["noCells"][1] += 2
        for name, size in zip(["xmx", "ymy", "zmz"], ["dx", "dy", "dz"]):
            dic[f"{name}_center"] = (dic[f"{name}"][1:] + dic[f"{name}"][:-1]) / 2.0
            dic[f"{size}"] = dic[f"{name}"][1:] - dic[f"{name}"][:-1]


def structured_handling_spe11a(dic):
    """
    Locate the geological positions in the tensor/cartesian grid for spe11a

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    sensor1, sensor2, centers, corners = [], [], [], []
    for k in range(dic["noCells"][2]):
        for i in range(dic["noCells"][0]):
            fgl = pd.Series(
                (
                    (dic["cxc1"] - dic["xmx_center"][i]) ** 2
                    + (dic["czc1"] - dic["zmz_center"][k]) ** 2
                )
            ).argmin()
            sensor1.append(
                (dic["xmx_center"][i] - dic["sensors"][0][0]) ** 2
                + (dic["zmz_center"][k] + dic["sensors"][0][2] - dic["dims"][2]) ** 2
            )
            sensor2.append(
                (dic["xmx_center"][i] - dic["sensors"][1][0]) ** 2
                + (dic["zmz_center"][k] + dic["sensors"][1][2] - dic["dims"][2]) ** 2
            )
            dic["satnum"].append(dic["ids_gmsh"][fgl][0])
            boxes(dic, dic["xmx_center"][i], dic["zmz_center"][k], i, dic["satnum"][-1])
            dic["permx"].append(f"{dic['rock'][int(dic['ids_gmsh'][fgl][0]) - 1][0]}")
            dic["poro"].append(f"{dic['rock'][int(dic['ids_gmsh'][fgl][0]) - 1][1]}")
            dic["disperc"].append(
                f"{dic['dispersion'][int(dic['ids_gmsh'][fgl][0])-1]}"
            )
            centers.append(
                f"{dic['xmx_center'][i]}, {dic['ymy_center'][0]}, {dic['zmz_center'][k]}"
            )
            corners.append(
                f"{dic['xmx'][i]}, {dic['dims'][2] -dic['zmz'][k]}, {dic['xmx'][i+1]}, "
                + f"{dic['dims'][2] -dic['zmz'][k]}, {dic['xmx'][i+1]}, "
                + f"{dic['dims'][2] -dic['zmz'][k+1]}, {dic['xmx'][i]}, "
                + f"{dic['dims'][2] -dic['zmz'][k+1]}"
            )
    dic["pop1"] = pd.Series(sensor1).argmin()
    dic["pop2"] = pd.Series(sensor2).argmin()
    dic["fipnum"][dic["pop1"]] = "8"
    dic["fipnum"][dic["pop2"]] = "9"
    sensors(dic)
    wells(dic)
    with open(
        f"{dic['exe']}/{dic['fol']}/deck/centers.txt",
        "w",
        encoding="utf8",
    ) as file:
        file.write("\n".join(centers))
    with open(
        f"{dic['exe']}/{dic['fol']}/deck/corners.txt",
        "w",
        encoding="utf8",
    ) as file:
        file.write("\n".join(corners))


def structured_handling_spe11bc(dic):
    """
    Locate the geological positions in the tensor/cartesian grid for the spe11b/c

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    sensor1, sensor2, centers, corners, pv_l = [], [], [], [], 0
    for k in range(dic["noCells"][2]):
        for i in range(dic["noCells"][0]):
            fgl = pd.Series(
                (
                    (dic["cxc1"] - dic["xmx_center"][i]) ** 2
                    + (dic["czc1"] - dic["zmz_center"][k]) ** 2
                )
            ).argmin()
            sensor1.append(
                (dic["xmx_center"][i] - dic["sensors"][0][0]) ** 2
                + (dic["ymy_center"][0] - dic["sensors"][0][1]) ** 2
                + (dic["zmz_center"][k] + dic["sensors"][0][2] - dic["dims"][2]) ** 2
            )
            sensor2.append(
                (dic["xmx_center"][i] - dic["sensors"][1][0]) ** 2
                + (dic["ymy_center"][0] - dic["sensors"][1][1]) ** 2
                + (dic["zmz_center"][k] + dic["sensors"][1][2] - dic["dims"][2]) ** 2
            )
            z_c = dic["zmz_center"][k]
            if dic["spe11"] == "spe11c":
                z_c -= map_z(dic, 0)
            dic["satnum"].append(dic["ids_gmsh"][fgl][0])
            boxes(dic, dic["xmx_center"][i], z_c, i, dic["satnum"][-1])
            dic["permx"].append(f"{dic['rock'][int(dic['ids_gmsh'][fgl][0]) - 1][0]}")
            poro = f"{dic['rock'][int(dic['ids_gmsh'][fgl][0]) - 1][1]}"
            dic["poro"].append(poro)
            pv = float(poro) * (dic["pvAdded"] + dic["widthBuffer"])
            dic["thconr"].append(f"{dic['rockCond'][int(dic['ids_gmsh'][fgl][0])-1]}")
            dic["disperc"].append(
                f"{dic['dispersion'][int(dic['ids_gmsh'][fgl][0])-1]}"
            )
            if i == 0 and (
                int(dic["ids_gmsh"][fgl][0]) != 1 and int(dic["ids_gmsh"][fgl][0]) != 7
            ):
                dic["porv"].append(
                    f"PORV {pv*dic['dy'][0]*dic['dz'][k]} 1 1 1 1 {k+1} {k+1} /"
                )
                pv_l = pv
            elif i == dic["noCells"][0] - 1 and (
                int(dic["ids_gmsh"][fgl][0]) != 1 and int(dic["ids_gmsh"][fgl][0]) != 7
            ):
                dic["porv"].append(
                    f"PORV {pv*dic['dy'][0]*dic['dz'][k]} {dic['noCells'][0]} "
                    + f"{dic['noCells'][0]} 1 1 {k+1} {k+1} /"
                )
            centers.append(
                f"{dic['xmx_center'][i]}, {dic['ymy_center'][0]}, {dic['zmz_center'][k]}"
            )
            corners.append(
                f"{dic['xmx'][i]}, {dic['dims'][2] -dic['zmz'][k]}, {dic['xmx'][i+1]}, "
                + f"{dic['dims'][2] -dic['zmz'][k]}, {dic['xmx'][i+1]}, "
                + f"{dic['dims'][2] -dic['zmz'][k+1]}, {dic['xmx'][i]}, "
                + f"{dic['dims'][2] -dic['zmz'][k+1]}"
            )
        for j in range(dic["noCells"][1] - 1):
            for names in ["satnum", "poro", "permx", "disperc", "thconr"]:
                dic[f"{names}"].extend(dic[f"{names}"][-dic["noCells"][0] :])
            for i_i in range(dic["noCells"][0]):
                sensor1.append(
                    (dic["xmx_center"][i_i] - dic["sensors"][0][0]) ** 2
                    + (dic["ymy_center"][j + 1] - dic["sensors"][0][1]) ** 2
                    + (dic["zmz_center"][k] + dic["sensors"][0][2] - dic["dims"][2])
                    ** 2
                )
                sensor2.append(
                    (dic["xmx_center"][i_i] - dic["sensors"][1][0]) ** 2
                    + (dic["ymy_center"][j + 1] - dic["sensors"][1][1]) ** 2
                    + (dic["zmz_center"][k] + dic["sensors"][1][2] - dic["dims"][2])
                    ** 2
                )
                z_c = dic["zmz_center"][k]
                if dic["spe11"] == "spe11c":
                    z_c -= map_z(dic, j + 1)
                boxes(
                    dic,
                    dic["xmx_center"][i_i],
                    z_c,
                    i_i,
                    dic["satnum"][-dic["noCells"][0] + i_i],
                )
                if i_i == 0 and (
                    int(dic["satnum"][-dic["noCells"][0] + i_i]) != 1
                    and int(dic["satnum"][-dic["noCells"][0] + i_i]) != 7
                ):
                    dic["porv"].append(
                        f"PORV {pv_l*dic['dy'][j+1]*dic['dz'][k]} 1 1 "
                        + f"{j+2} {j+2} {k+1} {k+1} /"
                    )
                elif i_i == dic["noCells"][0] - 1 and (
                    int(dic["satnum"][-dic["noCells"][0] + i_i]) != 1
                    and int(dic["satnum"][-dic["noCells"][0] + i_i]) != 7
                ):
                    dic["porv"].append(
                        f"PORV {pv*dic['dy'][j+1]*dic['dz'][k]} {dic['noCells'][0]} "
                        + f"{dic['noCells'][0]} {j+2} {j+2} {k+1} {k+1} /"
                    )
    if dic["spe11"] == "spe11c":
        add_pv_fipnum_front_back(dic)
    dic["pop1"] = pd.Series(sensor1).argmin()
    dic["pop2"] = pd.Series(sensor2).argmin()
    dic["fipnum"][dic["pop1"]] = "8"
    dic["fipnum"][dic["pop2"]] = "9"
    sensors(dic)
    wells(dic)
    with open(
        f"{dic['exe']}/{dic['fol']}/deck/centers.txt",
        "w",
        encoding="utf8",
    ) as file:
        file.write("\n".join(centers))
    with open(
        f"{dic['exe']}/{dic['fol']}/deck/corners.txt",
        "w",
        encoding="utf8",
    ) as file:
        file.write("\n".join(corners))


def add_pv_fipnum_front_back(dic):
    """
    Add the buffer pore volume and bc labels also on the front and back boundaries.

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    for k in range(dic["noCells"][2]):
        for i in range(dic["noCells"][0] - 2):
            ind = i + 1 + k * dic["noCells"][0] * dic["noCells"][1]
            if int(dic["satnum"][ind]) != 1 and int(dic["satnum"][ind]) != 7:
                pv = float(dic["poro"][ind]) * (dic["pvAdded"] + dic["widthBuffer"])
                if dic["grid"] == "corner-point":
                    ind_xz = i + 1 + k * dic["noCells"][0]
                    dic["porv"].append(
                        f"PORV {pv*dic['d_x'][i+1]*dic['d_z'][ind_xz]} {i+2} {i+2} "
                        + f"1 1 {k+1} {k+1} /"
                    )
                    dic["porv"].append(
                        f"PORV {pv*dic['d_x'][i+1]*dic['d_z'][ind_xz]} {i+2} {i+2} "
                        + f"{dic['noCells'][1]} {dic['noCells'][1]} {k+1} {k+1} /"
                    )
                else:
                    dic["porv"].append(
                        f"PORV {pv*dic['dx'][i+1]*dic['dz'][k]} {i+2} {i+2} "
                        + f"1 1 {k+1} {k+1} /"
                    )
                    dic["porv"].append(
                        f"PORV {pv*dic['dx'][i+1]*dic['dz'][k]} {i+2} {i+2} "
                        + f"{dic['noCells'][1]} {dic['noCells'][1]} {k+1} {k+1} /"
                    )
            set_back_front_fipnums(dic, ind)


def set_back_front_fipnums(dic, ind):
    """
    For the front and back boundaries in spe11c:\n
    Box A: Fipnum 13\n
    Satnum 1 and Box A: Fipnum 14\n
    Box B: Fipnum 15\n
    Satnum 1 and Box B: Fipnum 16\n
    Box C: Fipnum 17\n
    Satnum 1 and Box C: Fipnum 18\n

    Args:
        dic (dict): Global dictionary
        ind (int): ID index for the property

    Returns:
        dic (dict): Modified global dictionary

    """
    if int(dic["fipnum"][ind]) == 2:
        dic["fipnum"][ind] = "13"
        dic["fipnum"][ind + dic["noCells"][0] * (dic["noCells"][1] - 1)] = "13"
    elif int(dic["fipnum"][ind]) == 5:
        dic["fipnum"][ind] = "14"
        dic["fipnum"][ind + dic["noCells"][0] * (dic["noCells"][1] - 1)] = "14"
    elif int(dic["fipnum"][ind]) == 3:
        dic["fipnum"][ind] = "15"
        dic["fipnum"][ind + dic["noCells"][0] * (dic["noCells"][1] - 1)] = "15"
    elif int(dic["fipnum"][ind]) == 6:
        dic["fipnum"][ind] = "16"
        dic["fipnum"][ind + dic["noCells"][0] * (dic["noCells"][1] - 1)] = "16"
    elif int(dic["fipnum"][ind]) == 4:
        dic["fipnum"][ind] = "17"
        dic["fipnum"][ind + dic["noCells"][0] * (dic["noCells"][1] - 1)] = "17"
    elif int(dic["fipnum"][ind]) == 12:
        dic["fipnum"][ind] = "18"
        dic["fipnum"][ind + dic["noCells"][0] * (dic["noCells"][1] - 1)] = "18"
    elif int(dic["satnum"][ind]) == 1:
        dic["fipnum"][ind] = "10"
        dic["fipnum"][ind + dic["noCells"][0] * (dic["noCells"][1] - 1)] = "10"
    else:
        dic["fipnum"][ind] = "11"
        dic["fipnum"][ind + dic["noCells"][0] * (dic["noCells"][1] - 1)] = "11"


def corner_point_handling_spe11a(dic):
    """
    Locate the geological positions in the corner-point grid for the spe11a

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    well1, well2, sensor1, sensor2, centers, corners = [], [], [], [], [], []
    dic["wellijk"] = [[] for _ in range(len(dic["wellCoord"]))]
    for i in range(dic["no_cells"]):
        get_cell_info(dic, i)
        fgl = pd.Series(
            ((dic["cxc1"] - dic["xyz"][0]) ** 2 + (dic["czc1"] - dic["xyz"][2]) ** 2)
        ).argmin()
        well1.append(
            (dic["wellCoord"][0][0] - dic["xyz"][0]) ** 2
            + (dic["wellCoord"][0][2] - dic["xyz"][2]) ** 2
        )
        well2.append(
            (dic["wellCoord"][1][0] - dic["xyz"][0]) ** 2
            + (dic["wellCoord"][1][2] - dic["xyz"][2]) ** 2
        )
        sensor1.append(
            (dic["xyz"][0] - dic["sensors"][0][0]) ** 2
            + (dic["xyz"][2] + dic["sensors"][0][2] - dic["dims"][2]) ** 2
        )
        sensor2.append(
            (dic["xyz"][0] - dic["sensors"][1][0]) ** 2
            + (dic["xyz"][2] + dic["sensors"][1][2] - dic["dims"][2]) ** 2
        )
        dic["satnum"].append(dic["ids_gmsh"][fgl][0])
        boxes(dic, dic["xyz"][0], dic["xyz"][2], dic["ijk"][0], dic["satnum"][-1])
        dic["permx"].append(f"{dic['rock'][int(dic['ids_gmsh'][fgl][0]) - 1][0]}")
        dic["poro"].append(f"{dic['rock'][int(dic['ids_gmsh'][fgl][0]) - 1][1]}")
        dic["disperc"].append(f"{dic['dispersion'][int(dic['ids_gmsh'][fgl][0])-1]}")
        centers.append(f"{dic['xyz'][0]}, {dic['ymy_center'][0]}, {dic['xyz'][2]}")
        corners.append(dic["corns"])
    dic["pop1"] = pd.Series(sensor1).argmin()
    dic["pop2"] = pd.Series(sensor2).argmin()
    dic["fipnum"][dic["pop1"]] = "8"
    dic["fipnum"][dic["pop2"]] = "9"
    idwell1 = pd.Series(well1).argmin()
    idwell2 = pd.Series(well2).argmin()
    if dic["use"] == "opm":
        well1ijk = dic["gridf"].ijk_from_global_index(idwell1)
        well2ijk = dic["gridf"].ijk_from_global_index(idwell2)
        dic["sensorijk"][0] = dic["gridf"].ijk_from_global_index(dic["pop1"])
        dic["sensorijk"][1] = dic["gridf"].ijk_from_global_index(dic["pop2"])
    else:
        well1ijk = dic["gridf"].get_ijk(global_index=idwell1)
        well2ijk = dic["gridf"].get_ijk(global_index=idwell2)
        dic["sensorijk"][0] = dic["gridf"].get_ijk(global_index=dic["pop1"])
        dic["sensorijk"][1] = dic["gridf"].get_ijk(global_index=dic["pop2"])
    dic["wellijk"][0] = [well1ijk[0] + 1, 1, well1ijk[2] + 1]
    dic["wellijk"][1] = [well2ijk[0] + 1, 1, well2ijk[2] + 1]
    with open(
        f"{dic['exe']}/{dic['fol']}/deck/centers.txt",
        "w",
        encoding="utf8",
    ) as file:
        file.write("\n".join(centers))
    with open(
        f"{dic['exe']}/{dic['fol']}/deck/corners.txt",
        "w",
        encoding="utf8",
    ) as file:
        file.write("\n".join(corners))


def get_cell_info(dic, i):
    """
    Get the cell center coordinate and ijk from the simulation cell

    Args:
        dic (dict): Global dictionary\n
        i (int): Global cell index

    Returns:
        dic (dict): Modified global dictionary

    """
    if dic["use"] == "opm":
        dic["ijk"] = dic["gridf"].ijk_from_global_index(i)
        vxyz = dic["gridf"].xyz_from_ijk(dic["ijk"][0], dic["ijk"][1], dic["ijk"][2])
        dic["corns"] = (
            f"{vxyz[0][0]}, {dic['dims'][2] -vxyz[2][0]}, {vxyz[0][1]}, "
            + f"{dic['dims'][2] -vxyz[2][1]}, {vxyz[0][5]}, {dic['dims'][2] -vxyz[2][5]},"
            + f" {vxyz[0][4]}, {dic['dims'][2] - vxyz[2][4]}"
        )
        pxyz = Polygon(
            [
                [vxyz[0][0], vxyz[2][0]],
                [vxyz[0][1], vxyz[2][1]],
                [vxyz[0][5], vxyz[2][5]],
                [vxyz[0][4], vxyz[2][4]],
                [vxyz[0][0], vxyz[2][0]],
            ]
        ).centroid.wkt
        dic["xyz"] = list(float(j) for j in pxyz[7:-1].split(" "))
        dic["xyz"].insert(1, 0.0)
    else:
        dic["xyz"] = dic["gridf"].get_xyz(global_index=i)
        dic["ijk"] = dic["gridf"].get_ijk(global_index=i)
        if "vxyz" not in dic:
            dic["vxyz"] = dic["gridf"].export_corners(dic["gridf"].export_index())
        dic["corns"] = (
            f"{dic['vxyz'][i][0]}, {dic['dims'][2] - dic['vxyz'][i][2]}, {dic['vxyz'][i][3]}, "
            + f"{dic['dims'][2] - dic['vxyz'][i][5]}, {dic['vxyz'][i][15]}, "
            + f"{dic['dims'][2] - dic['vxyz'][i][17]}, "
            + f"{dic['vxyz'][i][12]}, {dic['dims'][2] - dic['vxyz'][i][14]}"
        )


def corner_point_handling_spe11bc(dic):
    """
    Locate the geological positions in the corner-point grid for the spe11b/c

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    well1, well2, sensor1, sensor2, xtemp, ztemp, centers, corners, pv_l = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        0,
    )
    dic["wellijk"] = [[] for _ in range(len(dic["wellCoord"]))]
    for i in range(dic["no_cells"]):
        get_cell_info(dic, i)
        xtemp.append(dic["xyz"][0])
        ztemp.append(dic["xyz"][2])
        fgl = pd.Series(
            ((dic["cxc1"] - dic["xyz"][0]) ** 2 + (dic["czc1"] - dic["xyz"][2]) ** 2)
        ).argmin()
        well1.append(
            (dic["wellCoord"][0][0] - dic["xyz"][0]) ** 2
            + (dic["wellCoord"][0][2] - dic["xyz"][2]) ** 2
        )
        well2.append(
            (dic["wellCoord"][1][0] - dic["xyz"][0]) ** 2
            + (dic["wellCoord"][1][2] - dic["xyz"][2]) ** 2
        )
        sensor1.append(
            (dic["xyz"][0] - dic["sensors"][0][0]) ** 2
            + (dic["xyz"][2] + dic["sensors"][0][2] - dic["dims"][2]) ** 2
        )
        sensor2.append(
            (dic["xyz"][0] - dic["sensors"][1][0]) ** 2
            + (dic["xyz"][2] + dic["sensors"][1][2] - dic["dims"][2]) ** 2
        )
        z_c = dic["xyz"][2]
        if dic["spe11"] == "spe11c":
            z_c -= map_z(dic, dic["ijk"][1])
        dic["satnum"].append(dic["ids_gmsh"][fgl][0])
        boxes(dic, dic["xyz"][0], z_c, dic["ijk"][0], dic["satnum"][-1])
        dic["permx"].append(f"{dic['rock'][int(dic['ids_gmsh'][fgl][0]) - 1][0]}")
        poro = f"{dic['rock'][int(dic['ids_gmsh'][fgl][0]) - 1][1]}"
        dic["poro"].append(poro)
        pv = float(poro) * (dic["pvAdded"] + dic["widthBuffer"])
        dic["disperc"].append(f"{dic['dispersion'][int(dic['ids_gmsh'][fgl][0])-1]}")
        dic["thconr"].append(f"{dic['rockCond'][int(dic['ids_gmsh'][fgl][0])-1]}")
        if dic["ijk"][0] == 0 and (
            int(dic["ids_gmsh"][fgl][0]) != 1 and int(dic["ids_gmsh"][fgl][0]) != 7
        ):
            dic["porv"].append(
                f"PORV { pv*dic['d_y'][0]*dic['d_z'][i]} 1 1 1 1 "
                + f"{dic['ijk'][2]+1} {dic['ijk'][2]+1} /"
            )
            pv_l = pv
        elif dic["ijk"][0] == dic["noCells"][0] - 1 and (
            int(dic["ids_gmsh"][fgl][0]) != 1 and int(dic["ids_gmsh"][fgl][0]) != 7
        ):
            dic["porv"].append(
                f"PORV {pv*dic['d_y'][0]*dic['d_z'][i]} {dic['noCells'][0]} "
                + f"{dic['noCells'][0]} 1 1 {dic['ijk'][2]+1} {dic['ijk'][2]+1} /"
            )
        centers.append(f"{dic['xyz'][0]}, {dic['ymy_center'][0]}, {dic['xyz'][2]}")
        corners.append(dic["corns"])
        if dic["ijk"][0] > 0 and dic["ijk"][0] == dic["noCells"][0] - 1:
            for j in range(dic["noCells"][1] - 1):
                for names in ["satnum", "poro", "permx", "disperc", "thconr"]:
                    dic[f"{names}"].extend(dic[f"{names}"][-dic["noCells"][0] :])
                for i_i in range(dic["noCells"][0]):
                    z_c = ztemp[i_i]
                    if dic["spe11"] == "spe11c":
                        z_c -= map_z(dic, j + 1)
                    boxes(
                        dic,
                        xtemp[i_i],
                        z_c,
                        i_i,
                        dic["satnum"][-dic["noCells"][0] + i_i],
                    )
                    if i_i == 0 and (
                        int(dic["satnum"][-dic["noCells"][0] + i_i]) != 1
                        and int(dic["satnum"][-dic["noCells"][0] + i_i]) != 7
                    ):
                        dic["d_zl"] = dic["d_z"][-dic["noCells"][0] + 1 + i]
                        dic["porv"].append(
                            "PORV "
                            + f"{pv_l*dic['d_y'][j+1]*dic['d_zl']} 1 1 "
                            + f"{j+2} {j+2} {dic['ijk'][2]+1} {dic['ijk'][2]+1} /"
                        )
                    elif i_i == dic["noCells"][0] - 1 and (
                        int(dic["satnum"][-dic["noCells"][0] + i_i]) != 1
                        and int(dic["satnum"][-dic["noCells"][0] + i_i]) != 7
                    ):
                        dic["porv"].append(
                            f"PORV {pv*dic['d_y'][j+1]*dic['d_z'][i]} "
                            + f"{dic['noCells'][0]} {dic['noCells'][0]} {j+2} {j+2} "
                            + f"{dic['ijk'][2]+1} {dic['ijk'][2]+1} /"
                        )
            xtemp, ztemp = [], []
    if dic["spe11"] == "spe11c":
        add_pv_fipnum_front_back(dic)
    dic["pop1"] = pd.Series(sensor1).argmin()
    dic["pop2"] = pd.Series(sensor2).argmin()
    dic["well1"] = pd.Series(well1).argmin()
    dic["well2"] = pd.Series(well2).argmin()
    locate_wells_sensors(dic)
    with open(
        f"{dic['exe']}/{dic['fol']}/deck/centers.txt",
        "w",
        encoding="utf8",
    ) as file:
        file.write("\n".join(centers))
    with open(
        f"{dic['exe']}/{dic['fol']}/deck/corners.txt",
        "w",
        encoding="utf8",
    ) as file:
        file.write("\n".join(corners))


def locate_wells_sensors(dic):
    """
    Find the wells/sources and sensors ijk positions

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    if dic["use"] == "opm":
        well1ijk = dic["gridf"].ijk_from_global_index(dic["well1"])
        well2ijk = dic["gridf"].ijk_from_global_index(dic["well2"])
        dic["sensorijk"][0] = list(dic["gridf"].ijk_from_global_index(dic["pop1"]))
        dic["sensorijk"][1] = list(dic["gridf"].ijk_from_global_index(dic["pop2"]))
    else:
        well1ijk = dic["gridf"].get_ijk(global_index=dic["well1"])
        well2ijk = dic["gridf"].get_ijk(global_index=dic["well2"])
        dic["sensorijk"][0] = list(dic["gridf"].get_ijk(global_index=dic["pop1"]))
        dic["sensorijk"][1] = list(dic["gridf"].get_ijk(global_index=dic["pop2"]))
    dic["wellijk"][0] = [well1ijk[0] + 1, 1, well1ijk[2] + 1]
    dic["wellijk"][1] = [well2ijk[0] + 1, 1, well2ijk[2] + 1]
    if dic["spe11"] == "spe11c":
        dic["wellijkf"] = [[] for _ in range(len(dic["wellCoord"]))]
        dic["wellijkf"][0] = [dic["wellijk"][0][0], 1, dic["wellijk"][0][2]]
        dic["wellijkf"][1] = [dic["wellijk"][1][0], 1, dic["wellijk"][1][2]]
        dic["ymy_center"] = (np.array(dic["ymy"][1:]) + np.array(dic["ymy"][:-1])) / 2.0
        wji = pd.Series(abs(dic["wellCoord"][0][1] - dic["ymy_center"])).argmin() + 1
        wjf = pd.Series(abs(dic["wellCoordF"][0][1] - dic["ymy_center"])).argmin() + 1
        sj1 = pd.Series(abs(dic["sensors"][0][1] - dic["ymy_center"])).argmin()
        sj2 = pd.Series(abs(dic["sensors"][1][1] - dic["ymy_center"])).argmin()
        dic["sensorijk"][0][1] = sj1
        dic["sensorijk"][1][1] = sj2
        dic["wellijk"][0][1] = wji
        dic["wellijk"][1][1] = wji
        dic["wellijkf"][0][1] = wjf
        dic["wellijkf"][1][1] = wjf
        dic["wellkh"] = []
        z_centers = []
        for k in range(dic["noCells"][2]):
            get_cell_info(dic, well1ijk[0] + k * dic["noCells"][0])
            z_centers.append(dic["xyz"][2])
        for j in range(dic["wellijk"][0][1], dic["wellijkf"][0][1] + 1):
            midpoints = z_centers - map_z(dic, j - 1)
            dic["wellkh"].append(
                pd.Series(
                    abs(
                        dic["wellCoord"][0][2]
                        - map_z(dic, dic["wellijk"][0][1] - 1)
                        - midpoints
                    )
                ).argmin()
                + 1
            )
    dic["fipnum"][
        dic["sensorijk"][0][0]
        + dic["sensorijk"][0][1] * dic["noCells"][0]
        + dic["sensorijk"][0][2] * dic["noCells"][0] * dic["noCells"][1]
    ] = "8"
    dic["fipnum"][
        dic["sensorijk"][1][0]
        + dic["sensorijk"][1][1] * dic["noCells"][0]
        + dic["sensorijk"][1][2] * dic["noCells"][0] * dic["noCells"][1]
    ] = "9"


def boxes(dic, x_c, z_c, idx, satnum):
    """
    Find the global indices for the different boxes for the report data

    Args:
        dic (dict): Global dictionary\n
        x_c (float): x-position of the cell center\n
        z_c (float): z-position of the cell center\n
        idx (int): i index of the cell position\n
        satnum (int): Number of the facie in the cell

    Returns:
        dic (dict): Modified global dictionary

    """
    if (
        (dic["dims"][2] + dic["maxelevation"] - z_c >= dic["boxb"][0][2])
        & (dic["dims"][2] + dic["maxelevation"] - z_c <= dic["boxb"][1][2])
        & (x_c >= dic["boxb"][0][0])
        & (x_c <= dic["boxb"][1][0])
    ):
        check_facie1(dic, satnum, "6", "3")
    elif (
        (dic["dims"][2] + dic["maxelevation"] - z_c >= dic["boxc"][0][2])
        & (dic["dims"][2] + dic["maxelevation"] - z_c <= dic["boxc"][1][2])
        & (x_c >= dic["boxc"][0][0])
        & (x_c <= dic["boxc"][1][0])
    ):
        check_facie1(dic, satnum, "12", "4")
    elif (
        (dic["dims"][2] + dic["maxelevation"] - z_c >= dic["boxa"][0][2])
        & (dic["dims"][2] + dic["maxelevation"] - z_c <= dic["boxa"][1][2])
        & (x_c >= dic["boxa"][0][0])
        & (x_c <= dic["boxa"][1][0])
    ):
        check_facie1(dic, satnum, "5", "2")
    elif dic["spe11"] != "spe11a" and idx in (0, dic["noCells"][0] - 1):
        check_facie1(dic, satnum, "10", "11")
    elif satnum == "1":
        dic["fipnum"].append("7")
    else:
        dic["fipnum"].append("1")


def check_facie1(dic, satnum, numa, numb):
    """
    Handle the overlaping with facie 1

    Args:
        dic (dict): Global dictionary\n
        satnum (int): Number of the facie in the cell\n
        numa (int): Fipnum to assign to the cell if it overlaps with Facie 1\n
        numb (int): Fipnum to assign to the cell otherwise.

    Returns:
        dic (dict): Modified global dictionary

    """
    if satnum == "1":
        dic["fipnum"].append(numa)
    else:
        dic["fipnum"].append(numb)


def positions(dic):
    """
    Function to locate sand and well positions

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    dic["sensorijk"] = [[] for _ in range(len(dic["sensors"]))]
    getfacies(dic)
    for names in ["satnum", "poro", "permx", "thconr", "fipnum", "disperc", "porv"]:
        dic[f"{names}"] = []
    if dic["grid"] == "corner-point":
        if dic["use"] == "opm":
            dic["gridf"] = OpmGrid(f"{dic['exe']}/{dic['fol']}/flow/INITIAL.EGRID")
            dic["initf"] = OpmFile(f"{dic['exe']}/{dic['fol']}/flow/INITIAL.INIT")
            dic["no_cells"] = (
                dic["gridf"].dimension[0]
                * dic["gridf"].dimension[1]
                * dic["gridf"].dimension[2]
            )
            dic["pv"] = dic["initf"]["PORV"]
            dic["actind"] = list(i for i in range(dic["no_cells"]) if dic["pv"][i] > 0)
            dic["d_z"] = np.array([0.0] * dic["no_cells"])
            dic["d_z"][dic["actind"]] = list(dic["initf"]["DZ"])
        else:
            dic["gridf"] = Grid(f"{dic['exe']}/{dic['fol']}/flow/INITIAL.EGRID")
            dic["initf"] = ResdataFile(f"{dic['exe']}/{dic['fol']}/flow/INITIAL.INIT")
            dic["actnum"] = list(dic["gridf"].export_actnum())
            dic["no_cells"] = dic["gridf"].nx * dic["gridf"].ny * dic["gridf"].nz
            dic["actind"] = list(i for i, act in enumerate(dic["actnum"]) if act == 1)
            dic["d_z"] = np.array([0.0] * dic["no_cells"])
            dic["d_z"][dic["actind"]] = list(dic["initf"].iget_kw("DZ")[0])
        if dic["spe11"] == "spe11a":
            corner_point_handling_spe11a(dic)
        else:
            corner_point_handling_spe11bc(dic)
    else:
        if dic["spe11"] == "spe11a":
            structured_handling_spe11a(dic)
        else:
            structured_handling_spe11bc(dic)
    np.savetxt(
        f"{dic['exe']}/{dic['fol']}/deck/ycenters.txt", dic["ymy_center"], fmt="%.8E"
    )


def sensors(dic):
    """
    Find the i,j,k sensor indices

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    for j, _ in enumerate(dic["sensors"]):
        for sensor_coord, axis in zip(dic["sensors"][j], ["xmx", "ymy", "zmz"]):
            if axis == "zmz":
                dic["sensorijk"][j].append(
                    pd.Series(
                        abs(dic["dims"][2] - sensor_coord - dic[f"{axis}_center"])
                    ).argmin()
                )
            else:
                dic["sensorijk"][j].append(
                    pd.Series(abs(sensor_coord - dic[f"{axis}_center"])).argmin()
                )


def wells(dic):
    """
    Function to find the wells/sources index

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    dic["wellijk"] = [[] for _ in range(len(dic["wellCoord"]))]
    if dic["spe11"] != "spe11c":
        for j, _ in enumerate(dic["wellCoord"]):
            for well_coord, axis in zip(dic["wellCoord"][j], ["xmx", "ymy", "zmz"]):
                dic["wellijk"][j].append(
                    pd.Series(abs(well_coord - dic[f"{axis}_center"])).argmin() + 1
                )
    else:
        dic["wellijkf"] = [[] for _ in range(len(dic["wellCoord"]))]
        for j, _ in enumerate(dic["wellCoord"]):
            for k, (well_coord, axis) in enumerate(
                zip(dic["wellCoord"][j][:2], ["xmx", "ymy"])
            ):
                dic["wellijk"][j].append(
                    pd.Series(abs(well_coord - dic[f"{axis}_center"])).argmin() + 1
                )
                dic["wellijkf"][j].append(
                    pd.Series(
                        abs(dic["wellCoordF"][j][k] - dic[f"{axis}_center"])
                    ).argmin()
                    + 1
                )
            if j == 0:
                well_coord = dic["wellCoord"][j][2]
                midpoints = dic["zmz_center"] - map_z(dic, dic["wellijk"][j][1] - 1)
                dic["wellijk"][j].append(
                    pd.Series(abs(well_coord - midpoints)).argmin() + 1
                )
            else:
                well_coord = dic["wellCoord"][j][2]
                midpoints = dic["zmz_center"]
                dic["wellijk"][j].append(
                    pd.Series(abs(well_coord - midpoints)).argmin() + 1
                )
    if dic["spe11"] == "spe11c":
        dic["wellkh"] = []
        for j in range(dic["wellijk"][0][1], dic["wellijkf"][0][1] + 1):
            midpoints = dic["zmz_center"] - map_z(dic, j - 1)
            dic["wellkh"].append(
                pd.Series(
                    abs(
                        dic["wellCoord"][0][2]
                        - map_z(dic, dic["wellijk"][0][1] - 1)
                        - midpoints
                    )
                ).argmin()
                + 1
            )


def map_z(dic, j):
    """
    Function to return the z position of the parabola for the wells

    Args:
        dic (dict): Global dictionary\n
        j : Cell id along the y axis

    Returns:
        z: Position
        dic (dict): Global dictionary

    """
    z_pos = (
        -dic["elevation"]
        + dic["elevation"]
        * (1.0 - (dic["ymy_center"][j] / (0.5 * dic["dims"][1]) - 1) ** 2.0)
        - dic["ymy_center"][j] * dic["backElevation"] / dic["dims"][1]
    )
    return z_pos


def getfacies(dic):
    """
    Function to read the reference Gmsh file

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    for ent in ["cxc1", "czc1", "ids_gmsh"]:
        dic[ent] = []
    # Read the points for the reference mesh file
    lol, points, simplex, centr, centrxz = ([], [], [], [], [])
    h_f_c, l_f_c = 0, 0

    if dic["spe11"] == "spe11a":
        h_ref = 1.2
        l_ref = 2.8
    else:
        h_ref = 1200.0
        l_ref = 8400.0
    with open(
        f"{dic['pat']}/reference_mesh/facies_coordinates.msh",
        "r",
        encoding="utf8",
    ) as file:
        for row in csv.reader(file, delimiter=" "):
            lol.append(row)
    for i in range(int(lol[14][0])):
        points.append(lol[15 + i][1:3])
        l_f_c, h_f_c = max(l_f_c, float(lol[15 + i][1])), max(
            h_f_c, float(lol[15 + i][2])
        )

    for i in range(int(lol[int(lol[14][0]) + 17][0])):
        dic["ids_gmsh"].append(lol[int(lol[14][0]) + 18 + i][3])
        simplex.append(lol[int(lol[14][0]) + 18 + i][5:])
        simplex[-1].insert(0, int(lol[int(lol[14][0]) + 18 + i][1]))

    for i, _ in enumerate(simplex):
        if simplex[i][0] == 2:
            centr.append(
                Polygon(
                    [
                        [
                            float(points[int(simplex[i][1]) - 1][0]) * l_ref / l_f_c,
                            dic["dims"][2]
                            - float(points[int(simplex[i][1]) - 1][1]) * h_ref / h_f_c,
                        ],
                        [
                            float(points[int(simplex[i][2]) - 1][0]) * l_ref / l_f_c,
                            dic["dims"][2]
                            - float(points[int(simplex[i][2]) - 1][1]) * h_ref / h_f_c,
                        ],
                        [
                            float(points[int(simplex[i][3]) - 1][0]) * l_ref / l_f_c,
                            dic["dims"][2]
                            - float(points[int(simplex[i][3]) - 1][1]) * h_ref / h_f_c,
                        ],
                        [
                            float(points[int(simplex[i][1]) - 1][0]) * l_ref / l_f_c,
                            dic["dims"][2]
                            - float(points[int(simplex[i][1]) - 1][1]) * h_ref / h_f_c,
                        ],
                    ]
                ).centroid.wkt
            )
        else:
            centr.append(
                Polygon(
                    [
                        [
                            float(points[int(simplex[i][1]) - 1][0]) * l_ref / l_f_c,
                            dic["dims"][2]
                            - float(points[int(simplex[i][1]) - 1][1]) * h_ref / h_f_c,
                        ],
                        [
                            float(points[int(simplex[i][2]) - 1][0]) * l_ref / l_f_c,
                            dic["dims"][2]
                            - float(points[int(simplex[i][2]) - 1][1]) * h_ref / h_f_c,
                        ],
                        [
                            float(points[int(simplex[i][3]) - 1][0]) * l_ref / l_f_c,
                            dic["dims"][2]
                            - float(points[int(simplex[i][3]) - 1][1]) * h_ref / h_f_c,
                        ],
                        [
                            float(points[int(simplex[i][4]) - 1][0]) * l_ref / l_f_c,
                            dic["dims"][2]
                            - float(points[int(simplex[i][4]) - 1][1]) * h_ref / h_f_c,
                        ],
                        [
                            float(points[int(simplex[i][1]) - 1][0]) * l_ref / l_f_c,
                            dic["dims"][2]
                            - float(points[int(simplex[i][1]) - 1][1]) * h_ref / h_f_c,
                        ],
                    ]
                ).centroid.wkt
            )
        centrxz.append([float(j) for j in centr[i][7:-1].split(" ")])
        dic["cxc1"].append(centrxz[i][0])
        dic["czc1"].append(centrxz[i][1])
    dic["cxc1"] = np.array(dic["cxc1"])
    dic["czc1"] = np.array(dic["czc1"])


def get_lines(dic):
    """
    Read the points in the z-surface lines

     Args:
        dic (dict): Global dictionary

     Returns:
        lines (list): List with the coordinates of the geological lines

    """
    with open(
        f"{dic['pat']}/reference_mesh/lines_coordinates.geo",
        "r",
        encoding="utf8",
    ) as file:
        lol = file.readlines()
    lines = []
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
                            float(points[0]) * dic["dims"][0] / 2.8,
                            (1.2 - float(points[1]) - float(points[2][:-3]))
                            * dic["dims"][2]
                            / 1.2,
                        ]
                    )
                    break
        else:
            newline = True
    return lines


def corner(dic):
    """
    Create a SPE11 corner-point grid

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    # Read the points for the .geo gmsh file
    lines = get_lines(dic)
    dic["xcor"], dic["zcor"] = [], []
    dic["xmx"] = [0.0]
    for i, n_x in enumerate(dic["x_n"]):
        for j in range(n_x):
            dic["xmx"].append((i + (j + 1.0) / n_x) * dic["dims"][0] / len(dic["x_n"]))
    dic["ymy"] = [0.0]
    for i, n_y in enumerate(dic["y_n"]):
        for j in range(n_y):
            dic["ymy"].append((i + (j + 1.0) / n_y) * dic["dims"][1] / len(dic["y_n"]))
    if (dic["spe11"] == "spe11b" or dic["spe11"] == "spe11c") and 1.1 * dic[
        "widthBuffer"
    ] < dic["xmx"][1]:
        dic["xmx"] = np.insert(dic["xmx"], 1, dic["widthBuffer"])
        dic["xmx"] = np.insert(
            dic["xmx"], len(dic["xmx"]) - 1, dic["xmx"][-1] - dic["widthBuffer"]
        )
    if dic["spe11"] == "spe11c" and 1.1 * dic["widthBuffer"] < dic["ymy"][1]:
        dic["ymy"] = np.insert(dic["ymy"], 1, dic["widthBuffer"])
        dic["ymy"] = np.insert(
            dic["ymy"], len(dic["ymy"]) - 1, dic["ymy"][-1] - dic["widthBuffer"]
        )
    dic["noCells"][1] = len(dic["ymy"]) - 1
    for xcor in dic["xmx"]:
        for _, lcor in enumerate(lines):
            dic["xcor"].append(xcor)
            idx = pd.Series([abs(ii[0] - xcor) for ii in lcor]).argmin()
            if lcor[idx][0] < xcor:
                dic["zcor"].append(
                    lcor[idx][1]
                    + (
                        (lcor[idx + 1][1] - lcor[idx][1])
                        / (lcor[idx + 1][0] - lcor[idx][0])
                    )
                    * (xcor - lcor[idx][0])
                )
            else:
                dic["zcor"].append(
                    lcor[idx - 1][1]
                    + (
                        (lcor[idx][1] - lcor[idx - 1][1])
                        / (lcor[idx][0] - lcor[idx - 1][0])
                    )
                    * (xcor - lcor[idx - 1][0])
                )
    res = list(filter(lambda i: i == dic["zcor"][-1], dic["zcor"]))[0]
    n_z = dic["zcor"].index(res)
    res = list(filter(lambda i: i > 0, dic["xcor"]))[0]
    n_x = round(len(dic["xcor"]) / dic["xcor"].index(res)) - 1
    dic["noCells"][0], dic["noCells"][2] = n_x, n_z
    # Refine the grid
    dic["xcor"], dic["zcor"], dic["noCells"][0], dic["noCells"][2] = refinement_z(
        dic["xcor"], dic["zcor"], dic["noCells"][0], dic["noCells"][2], dic["z_n"]
    )
    dic["xmx"] = np.array(dic["xmx"])
    dic["ymy_center"] = 0.5 * (np.array(dic["ymy"])[1:] + np.array(dic["ymy"])[:-1])
    dic["d_y"] = np.array(dic["ymy"])[1:] - np.array(dic["ymy"])[:-1]
    dic["d_x"] = np.array(dic["xmx"])[1:] - np.array(dic["xmx"])[:-1]


def refinement_z(xci, zci, ncx, ncz, znr):
    """
    Refinment of the grid in the z-dir

    Args:
        xci (list): Floats with the x-coordinates of the cell corners\n
        zci (list): Floats with the z-coordinates of the cell corners\n
        ncx (int): Number of cells in the x-dir\n
        ncz (int): Number of cells in the z-dir\n
        znr (list): Integers with the number of z-refinments per cell

    Returns:
        xcr (list): Floats with the new x-coordinates of the cell corners\n
        zcr (list): Floats with the new z-coordinates of the cell corners\n
        ncx (int): New number of cells in the x-dir\n
        ncz (int): New number of cells in the z-dir

    """
    xcr, zcr = [], []
    for j in range(ncx + 1):
        zcr.append(zci[j * (ncz + 1)])
        xcr.append(xci[j * (ncz + 1)])
        for i in range(ncz):
            for k in range(znr[i]):
                alp = np.arange(1.0 / znr[i], 1.0 + 1.0 / znr[i], 1.0 / znr[i]).tolist()
                zcr.append(
                    zci[j * (ncz + 1) + i]
                    + (zci[j * (ncz + 1) + i + 1] - zci[j * (ncz + 1) + i]) * alp[k]
                )
                xcr.append(
                    xci[j * (ncz + 1) + i]
                    + (xci[j * (ncz + 1) + i + 1] - xci[j * (ncz + 1) + i]) * alp[k]
                )
    res = list(filter(lambda i: i > 0, xcr))[0]
    ncx = round(len(xcr) / xcr.index(res)) - 1
    res = list(filter(lambda i: i == zcr[-1], zcr))[0]
    ncz = zcr.index(res)
    return xcr, zcr, ncx, ncz
