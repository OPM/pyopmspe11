# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT
# pylint: disable=C0302, R0912, R0914, R0915, E1102

"""
Utiliy function for the grid and locations in the geological models.
"""

import csv
import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon
from alive_progress import alive_bar
from pyopmspe11.utils.writefile import create_corner_point_grid


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
    sensor1, sensor2 = [], []
    with alive_bar(dic["noCells"][0] * dic["noCells"][2]) as bar_animation:
        for k in range(dic["noCells"][2]):
            for i in range(dic["noCells"][0]):
                bar_animation()
                n = -1
                if dic["zmz_center"][k] > dic["zmidbot"]:
                    order = dic["under_order"]
                elif dic["zmz_center"][k] > dic["ztopbot"]:
                    order = dic["bottom_order"]
                else:
                    order = dic["top_order"]
                for ind in order:
                    if dic["polygons"][ind].contains(
                        Point(dic["xmx_center"][i], dic["zmz_center"][k])
                    ):
                        n = dic["facies"][ind]
                        break
                sensor1.append(
                    (dic["xmx_center"][i] - dic["sensors"][0][0]) ** 2
                    + (dic["zmz_center"][k] + dic["sensors"][0][2] - dic["dims"][2])
                    ** 2
                )
                sensor2.append(
                    (dic["xmx_center"][i] - dic["sensors"][1][0]) ** 2
                    + (dic["zmz_center"][k] + dic["sensors"][1][2] - dic["dims"][2])
                    ** 2
                )
                dic["fluxnum"].append(str(n))
                boxes(
                    dic,
                    dic["xmx_center"][i],
                    dic["zmz_center"][k],
                    i,
                    dic["fluxnum"][-1],
                )
    dic["pop1"] = pd.Series(sensor1).argmin()
    dic["pop2"] = pd.Series(sensor2).argmin()
    dic["fipnum"][dic["pop1"]] = "8"
    dic["fipnum"][dic["pop2"]] = "9"
    sensors(dic)
    wells(dic)


def structured_handling_spe11bc(dic):
    """
    Locate the geological positions in the tensor/cartesian grid for the spe11b/c

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    sensor1, sensor2, pv_l = [], [], 0
    with alive_bar(dic["noCells"][0] * dic["noCells"][2]) as bar_animation:
        for k in range(dic["noCells"][2]):
            for i in range(dic["noCells"][0]):
                bar_animation()
                n = -1
                if dic["zmz_center"][k] > dic["zmidbot"]:
                    order = dic["under_order"]
                elif dic["zmz_center"][k] > dic["ztopbot"]:
                    order = dic["bottom_order"]
                else:
                    order = dic["top_order"]
                for ind in order:
                    if dic["polygons"][ind].contains(
                        Point(dic["xmx_center"][i], dic["zmz_center"][k])
                    ):
                        n = dic["facies"][ind]
                        break
                sensor1.append(
                    (dic["xmx_center"][i] - dic["sensors"][0][0]) ** 2
                    + (dic["ymy_center"][0] - dic["sensors"][0][1]) ** 2
                    + (dic["zmz_center"][k] + dic["sensors"][0][2] - dic["dims"][2])
                    ** 2
                )
                sensor2.append(
                    (dic["xmx_center"][i] - dic["sensors"][1][0]) ** 2
                    + (dic["ymy_center"][0] - dic["sensors"][1][1]) ** 2
                    + (dic["zmz_center"][k] + dic["sensors"][1][2] - dic["dims"][2])
                    ** 2
                )
                z_c = dic["zmz_center"][k]
                if dic["spe11"] == "spe11c":
                    z_c -= map_z(dic, 0)
                dic["fluxnum"].append(str(n))
                boxes(dic, dic["xmx_center"][i], z_c, i, dic["fluxnum"][-1])
                pv = dic["rock"][n - 1][1] * (dic["pvAdded"] + dic["widthBuffer"])
                if i == 0 and n not in (1, 7):
                    dic["porv"].append(
                        f"PORV {pv*dic['dy'][0]*dic['dz'][k]} 1 1 1 1 {k+1} {k+1} /"
                    )
                    pv_l = pv
                elif i == dic["noCells"][0] - 1 and n not in (1, 7):
                    dic["porv"].append(
                        f"PORV {pv*dic['dy'][0]*dic['dz'][k]} {dic['noCells'][0]} "
                        + f"{dic['noCells'][0]} 1 1 {k+1} {k+1} /"
                    )
            for j in range(dic["noCells"][1] - 1):
                dic["fluxnum"].extend(dic["fluxnum"][-dic["noCells"][0] :])
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
                        dic["fluxnum"][-dic["noCells"][0] + i_i],
                    )
                    if i_i == 0 and (
                        int(dic["fluxnum"][-dic["noCells"][0] + i_i]) != 1
                        and int(dic["fluxnum"][-dic["noCells"][0] + i_i]) != 7
                    ):
                        dic["porv"].append(
                            f"PORV {pv_l*dic['dy'][j+1]*dic['dz'][k]} 1 1 "
                            + f"{j+2} {j+2} {k+1} {k+1} /"
                        )
                    elif i_i == dic["noCells"][0] - 1 and (
                        int(dic["fluxnum"][-dic["noCells"][0] + i_i]) != 1
                        and int(dic["fluxnum"][-dic["noCells"][0] + i_i]) != 7
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
            if int(dic["fluxnum"][ind]) != 1 and int(dic["fluxnum"][ind]) != 7:
                pv = dic["rock"][int(dic["fluxnum"][ind]) - 1][1] * (
                    dic["pvAdded"] + dic["widthBuffer"]
                )
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
    Facie 1 and Box A: Fipnum 14\n
    Box B: Fipnum 15\n
    Facie 1 and Box B: Fipnum 16\n
    Box C: Fipnum 17\n
    Facie 1 and Box C: Fipnum 18\n

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
    elif int(dic["fluxnum"][ind]) == 1:
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
    well1, well2, sensor1, sensor2 = [], [], [], []
    dic["wellijk"] = [[] for _ in range(len(dic["wellCoord"]))]
    with alive_bar(dic["no_cells"]) as bar_animation:
        for i in range(dic["no_cells"]):
            bar_animation()
            i_x = int(i % dic["noCells"][0])
            k_z = int(np.floor(i / dic["noCells"][0]))
            n = -1
            if dic["xyz"][i][2] > dic["zmidbot"]:
                order = dic["under_order"]
            elif dic["xyz"][i][2] > dic["ztopbot"]:
                order = dic["bottom_order"]
            else:
                order = dic["top_order"]
            for ind in order:
                if dic["polygons"][ind].contains(
                    Point(dic["xyz"][i][0], dic["xyz"][i][2])
                ):
                    n = dic["facies"][ind]
                    break
            well1.append(
                (dic["wellCoord"][0][0] - dic["xyz"][i][0]) ** 2
                + (dic["wellCoord"][0][2] - dic["xyz"][i][2]) ** 2
            )
            well2.append(
                (dic["wellCoord"][1][0] - dic["xyz"][i][0]) ** 2
                + (dic["wellCoord"][1][2] - dic["xyz"][i][2]) ** 2
            )
            sensor1.append(
                (dic["xyz"][i][0] - dic["sensors"][0][0]) ** 2
                + (dic["xyz"][i][2] + dic["sensors"][0][2] - dic["dims"][2]) ** 2
            )
            sensor2.append(
                (dic["xyz"][i][0] - dic["sensors"][1][0]) ** 2
                + (dic["xyz"][i][2] + dic["sensors"][1][2] - dic["dims"][2]) ** 2
            )
            dic["fluxnum"].append(str(n))
            boxes(dic, dic["xyz"][i][0], dic["xyz"][i][2], i_x, dic["fluxnum"][-1])
    dic["pop1"] = pd.Series(sensor1).argmin()
    dic["pop2"] = pd.Series(sensor2).argmin()
    dic["fipnum"][dic["pop1"]] = "8"
    dic["fipnum"][dic["pop2"]] = "9"
    idwell1 = pd.Series(well1).argmin()
    idwell2 = pd.Series(well2).argmin()
    i_x = int(idwell1 % dic["noCells"][0])
    k_z = int(np.floor(idwell1 / dic["noCells"][0]))
    well1ijk = [i_x, 0, k_z]
    i_x = int(idwell2 % dic["noCells"][0])
    k_z = int(np.floor(idwell2 / dic["noCells"][0]))
    well2ijk = [i_x, 0, k_z]
    i_x = int(dic["pop1"] % dic["noCells"][0])
    k_z = int(np.floor(dic["pop1"] / dic["noCells"][0]))
    dic["sensorijk"][0] = [i_x, 0, k_z]
    i_x = int(dic["pop2"] % dic["noCells"][0])
    k_z = int(np.floor(dic["pop2"] / dic["noCells"][0]))
    dic["sensorijk"][1] = [i_x, 0, k_z]
    dic["wellijk"][0] = [well1ijk[0] + 1, 1, well1ijk[2] + 1]
    dic["wellijk"][1] = [well2ijk[0] + 1, 1, well2ijk[2] + 1]


def corner_point_handling_spe11bc(dic):
    """
    Locate the geological positions in the corner-point grid for the spe11b/c

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    well1, well2, sensor1, sensor2, xtemp, ztemp, pv_l = (
        [],
        [],
        [],
        [],
        [],
        [],
        0,
    )
    dic["wellijk"] = [[] for _ in range(len(dic["wellCoord"]))]
    with alive_bar(dic["no_cells"]) as bar_animation:
        for i in range(dic["no_cells"]):
            bar_animation()
            i_x = int(i % dic["noCells"][0])
            k_z = int(np.floor(i / dic["noCells"][0]))
            xtemp.append(dic["xyz"][i][0])
            ztemp.append(dic["xyz"][i][2])
            n = -1
            if dic["xyz"][i][2] > dic["zmidbot"]:
                order = dic["under_order"]
            elif dic["xyz"][i][2] > dic["ztopbot"]:
                order = dic["bottom_order"]
            else:
                order = dic["top_order"]
            for ind in order:
                if dic["polygons"][ind].contains(
                    Point(dic["xyz"][i][0], dic["xyz"][i][2])
                ):
                    n = dic["facies"][ind]
                    break
            well1.append(
                (dic["wellCoord"][0][0] - dic["xyz"][i][0]) ** 2
                + (dic["wellCoord"][0][2] - dic["xyz"][i][2]) ** 2
            )
            well2.append(
                (dic["wellCoord"][1][0] - dic["xyz"][i][0]) ** 2
                + (dic["wellCoord"][1][2] - dic["xyz"][i][2]) ** 2
            )
            sensor1.append(
                (dic["xyz"][i][0] - dic["sensors"][0][0]) ** 2
                + (dic["xyz"][i][2] + dic["sensors"][0][2] - dic["dims"][2]) ** 2
            )
            sensor2.append(
                (dic["xyz"][i][0] - dic["sensors"][1][0]) ** 2
                + (dic["xyz"][i][2] + dic["sensors"][1][2] - dic["dims"][2]) ** 2
            )
            z_c = dic["xyz"][i][2]
            if dic["spe11"] == "spe11c":
                z_c -= map_z(dic, 0)
            dic["fluxnum"].append(str(n))
            boxes(dic, dic["xyz"][i][0], z_c, i_x, dic["fluxnum"][-1])
            pv = dic["rock"][n - 1][1] * (dic["pvAdded"] + dic["widthBuffer"])
            if i_x == 0 and n not in (1, 7):
                dic["porv"].append(
                    f"PORV { pv*dic['d_y'][0]*dic['d_z'][i]} 1 1 1 1 "
                    + f"{k_z+1} {k_z+1} /"
                )
                pv_l = pv
            elif i_x == dic["noCells"][0] - 1 and n not in (1, 7):
                dic["porv"].append(
                    f"PORV {pv*dic['d_y'][0]*dic['d_z'][i]} {dic['noCells'][0]} "
                    + f"{dic['noCells'][0]} 1 1 {k_z+1} {k_z+1} /"
                )
            if i_x > 0 and i_x == dic["noCells"][0] - 1:
                for j in range(dic["noCells"][1] - 1):
                    dic["fluxnum"].extend(dic["fluxnum"][-dic["noCells"][0] :])
                    for i_i in range(dic["noCells"][0]):
                        z_c = ztemp[i_i]
                        if dic["spe11"] == "spe11c":
                            z_c -= map_z(dic, j + 1)
                        boxes(
                            dic,
                            xtemp[i_i],
                            z_c,
                            i_i,
                            dic["fluxnum"][-dic["noCells"][0] + i_i],
                        )
                        if i_i == 0 and (
                            int(dic["fluxnum"][-dic["noCells"][0] + i_i]) != 1
                            and int(dic["fluxnum"][-dic["noCells"][0] + i_i]) != 7
                        ):
                            dic["d_zl"] = dic["d_z"][-dic["noCells"][0] + 1 + i]
                            dic["porv"].append(
                                "PORV "
                                + f"{pv_l*dic['d_y'][j+1]*dic['d_zl']} 1 1 "
                                + f"{j+2} {j+2} {k_z+1} {k_z+1} /"
                            )
                        elif i_i == dic["noCells"][0] - 1 and (
                            int(dic["fluxnum"][-dic["noCells"][0] + i_i]) != 1
                            and int(dic["fluxnum"][-dic["noCells"][0] + i_i]) != 7
                        ):
                            dic["porv"].append(
                                f"PORV {pv*dic['d_y'][j+1]*dic['d_z'][i]} "
                                + f"{dic['noCells'][0]} {dic['noCells'][0]} {j+2} {j+2} "
                                + f"{k_z+1} {k_z+1} /"
                            )
                xtemp, ztemp = [], []
    if dic["spe11"] == "spe11c":
        add_pv_fipnum_front_back(dic)
    dic["pop1"] = pd.Series(sensor1).argmin()
    dic["pop2"] = pd.Series(sensor2).argmin()
    dic["well1"] = pd.Series(well1).argmin()
    dic["well2"] = pd.Series(well2).argmin()
    locate_wells_sensors(dic)


def locate_wells_sensors(dic):
    """
    Find the wells/sources and sensors ijk positions

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    i_x = int(dic["well1"] % dic["noCells"][0])
    k_z = int(np.floor(dic["well1"] / dic["noCells"][0]))
    well1ijk = [i_x, 0, k_z]
    i_x = int(dic["well2"] % dic["noCells"][0])
    k_z = int(np.floor(dic["well2"] / dic["noCells"][0]))
    well2ijk = [i_x, 0, k_z]
    i_x = int(dic["pop1"] % dic["noCells"][0])
    k_z = int(np.floor(dic["pop1"] / dic["noCells"][0]))
    dic["sensorijk"][0] = [i_x, 0, k_z]
    i_x = int(dic["pop2"] % dic["noCells"][0])
    k_z = int(np.floor(dic["pop2"] / dic["noCells"][0]))
    dic["sensorijk"][1] = [i_x, 0, k_z]
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
            z_centers.append(dic["xyz"][well1ijk[0] + k * dic["noCells"][0]][2])
        for j in range(dic["wellijk"][0][1], dic["wellijkf"][0][1] + 1):
            midpoints = z_centers - map_z(dic, j - 1) - dic["maxelevation"]
            dic["wellkh"].append(
                pd.Series(abs(dic["wellCoord"][0][2] - midpoints)).argmin() + 1
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


def boxes(dic, x_c, z_c, idx, fluxnum):
    """
    Find the global indices for the different boxes for the report data

    Args:
        dic (dict): Global dictionary\n
        x_c (float): x-position of the cell center\n
        z_c (float): z-position of the cell center\n
        idx (int): i index of the cell position\n
        fluxnum (int): Number of the facie in the cell

    Returns:
        dic (dict): Modified global dictionary

    """
    if (
        (dic["dims"][2] + dic["maxelevation"] - z_c >= dic["boxb"][0][2])
        & (dic["dims"][2] + dic["maxelevation"] - z_c <= dic["boxb"][1][2])
        & (x_c >= dic["boxb"][0][0])
        & (x_c <= dic["boxb"][1][0])
    ):
        check_facie1(dic, fluxnum, "6", "3")
    elif (
        (dic["dims"][2] + dic["maxelevation"] - z_c >= dic["boxc"][0][2])
        & (dic["dims"][2] + dic["maxelevation"] - z_c <= dic["boxc"][1][2])
        & (x_c >= dic["boxc"][0][0])
        & (x_c <= dic["boxc"][1][0])
    ):
        check_facie1(dic, fluxnum, "12", "4")
    elif (
        (dic["dims"][2] + dic["maxelevation"] - z_c >= dic["boxa"][0][2])
        & (dic["dims"][2] + dic["maxelevation"] - z_c <= dic["boxa"][1][2])
        & (x_c >= dic["boxa"][0][0])
        & (x_c <= dic["boxa"][1][0])
    ):
        check_facie1(dic, fluxnum, "5", "2")
    elif dic["spe11"] != "spe11a" and idx in (0, dic["noCells"][0] - 1):
        check_facie1(dic, fluxnum, "10", "11")
    elif fluxnum == "1":
        dic["fipnum"].append("7")
    else:
        dic["fipnum"].append("1")


def check_facie1(dic, fluxnum, numa, numb):
    """
    Handle the overlaping with facie 1

    Args:
        dic (dict): Global dictionary\n
        fluxnum (int): Number of the facie in the cell\n
        numa (int): Fipnum to assign to the cell if it overlaps with Facie 1\n
        numb (int): Fipnum to assign to the cell otherwise.

    Returns:
        dic (dict): Modified global dictionary

    """
    if fluxnum == "1":
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
    getpolygons(dic)
    for names in ["fluxnum", "fipnum", "porv"]:
        dic[f"{names}"] = []
    if dic["grid"] == "corner-point":
        if dic["spe11"] == "spe11a":
            corner_point_handling_spe11a(dic)
        else:
            corner_point_handling_spe11bc(dic)
    else:
        if dic["spe11"] == "spe11a":
            structured_handling_spe11a(dic)
        else:
            structured_handling_spe11bc(dic)


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
            midpoints = dic["zmz_center"] - map_z(dic, j - 1) - dic["maxelevation"]
            dic["wellkh"].append(
                pd.Series(abs(dic["wellCoord"][0][2] - midpoints)).argmin() + 1
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


def getpolygons(dic):
    """
    Function to create the polygons from the benchmark geo file

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    points = []
    lines = []
    curves = []
    dic["polygons"] = []
    dic["facies"] = []
    facie = 0
    if dic["spe11"] == "spe11a":
        h_ref = 1
        l_ref = 1
        dic["ztopbot"] = dic["dims"][2] - 0.644
        dic["zmidbot"] = dic["dims"][2] - 0.265
    else:
        h_ref = 1200.0 / 1.2
        l_ref = 8400.0 / 2.8
        dic["ztopbot"] = dic["dims"][2] - 0.644 * 1200.0 / 1.2
        dic["zmidbot"] = dic["dims"][2] - 0.265 * 1200.0 / 1.2
    with open(
        f"{dic['pat']}/reference_mesh/facies_coordinates.geo",
        "r",
        encoding="utf8",
    ) as file:
        for row in csv.reader(file, delimiter=" "):
            if row[0] == "//" and not points:
                continue
            if row[0][:5] == "Point":
                points.append(
                    [
                        l_ref * float(row[2][1:-1]),
                        dic["dims"][2] - h_ref * float(row[3][:-1]),
                    ]
                )
            elif row[0][:4] == "Line":
                lines.append([int(row[2][1:-1]), int(row[3][:-2])])
            elif row[0] == "Curve":
                dic["facies"].append(facie)
                tmp = []
                tmp.append(int(row[3][1:-1]))
                for col in row[4:-1]:
                    tmp.append(int(col[:-1]))
                tmp.append(int(row[-1][:-2]))
                curves.append(tmp)
            else:
                facie += 1
    for curve in curves:
        tmp = []
        tmp.append(
            [
                points[lines[curve[0] - 1][0] - 1][0],
                points[lines[curve[0] - 1][0] - 1][1],
            ]
        )
        tmp.append(
            [
                points[lines[curve[0] - 1][1] - 1][0],
                points[lines[curve[0] - 1][1] - 1][1],
            ]
        )
        for line in curve[1:]:
            if line < 0:
                tmp.append(
                    [
                        points[lines[abs(line) - 1][0] - 1][0],
                        points[lines[abs(line) - 1][0] - 1][1],
                    ]
                )
            else:
                tmp.append(
                    [
                        points[lines[line - 1][1] - 1][0],
                        points[lines[line - 1][1] - 1][1],
                    ]
                )
        tmp.append(tmp[0])
        dic["polygons"].append(Polygon(tmp))
    # Performance: Set the order to search in the polygons above and below
    dic["top_order"] = [
        0,
        4,
        5,
        6,
        9,
        10,
        11,
        15,
        16,
        17,
        18,
        22,
        23,
        24,
        27,
        28,
        30,
        7,
        13,
        12,
        20,
        1,
        8,
        19,
        21,
        14,
        2,
        3,
        25,
        26,
        29,
        31,
    ]
    dic["bottom_order"] = [
        12,
        20,
        1,
        8,
        19,
        21,
        14,
        2,
        3,
        25,
        26,
        29,
        31,
        0,
        4,
        5,
        6,
        9,
        10,
        11,
        15,
        16,
        17,
        18,
        22,
        23,
        24,
        27,
        28,
        30,
        7,
        13,
    ]
    dic["under_order"] = [
        25,
        26,
        31,
        29,
        12,
        20,
        1,
        8,
        19,
        21,
        14,
        2,
        3,
        0,
        4,
        5,
        6,
        9,
        10,
        11,
        15,
        16,
        17,
        18,
        22,
        23,
        24,
        27,
        28,
        30,
        7,
        13,
    ]


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
    xcoord, zcoord = [], []
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
            xcoord.append(xcor)
            idx = pd.Series([abs(ii[0] - xcor) for ii in lcor]).argmin()
            if lcor[idx][0] < xcor:
                zcoord.append(
                    lcor[idx][1]
                    + (
                        (lcor[idx + 1][1] - lcor[idx][1])
                        / (lcor[idx + 1][0] - lcor[idx][0])
                    )
                    * (xcor - lcor[idx][0])
                )
            else:
                zcoord.append(
                    lcor[idx - 1][1]
                    + (
                        (lcor[idx][1] - lcor[idx - 1][1])
                        / (lcor[idx][0] - lcor[idx - 1][0])
                    )
                    * (xcor - lcor[idx - 1][0])
                )
    res = list(filter(lambda i: i == zcoord[-1], zcoord))[0]
    n_z = zcoord.index(res)
    res = list(filter(lambda i: i > 0, xcoord))[0]
    n_x = round(len(xcoord) / xcoord.index(res)) - 1
    dic["noCells"][0], dic["noCells"][2] = n_x, n_z
    # Refine the grid
    xcoord, zcoord, dic["noCells"][0], dic["noCells"][2] = refinement_z(
        xcoord, zcoord, dic["noCells"][0], dic["noCells"][2], dic["z_n"]
    )
    dic["xmx"] = np.array(dic["xmx"])
    dic["ymy_center"] = 0.5 * (np.array(dic["ymy"])[1:] + np.array(dic["ymy"])[:-1])
    dic["d_y"] = np.array(dic["ymy"])[1:] - np.array(dic["ymy"])[:-1]
    dic["d_x"] = np.array(dic["xmx"])[1:] - np.array(dic["xmx"])[:-1]
    dic["no_cells"] = dic["noCells"][0] * dic["noCells"][2]
    create_corner_point_grid(dic, xcoord, zcoord)
    dic["xyz"] = []
    dic["d_z"] = []
    for k in range(dic["noCells"][2]):
        for i in range(dic["noCells"][0]):
            n = (i * (dic["noCells"][2] + 1)) + k
            m = ((i + 1) * (dic["noCells"][2] + 1)) + k
            poly = Polygon(
                [
                    [xcoord[n], zcoord[n]],
                    [xcoord[m], zcoord[m]],
                    [xcoord[m + 1], zcoord[m + 1]],
                    [xcoord[n + 1], zcoord[n + 1]],
                    [xcoord[n], zcoord[n]],
                ]
            )
            pxz = poly.centroid.wkt
            pxz = list(float(j) for j in pxz[7:-1].split(" "))
            dic["xyz"].append([pxz[0], 0, pxz[1]])
            dic["d_z"].append(poly.area / (xcoord[m] - xcoord[n]))


def refinement_z(xci, zci, ncx, ncz, znr):
    """
    Refinement of the grid in the z-dir

    Args:
        xci (list): Floats with the x-coordinates of the cell corners\n
        zci (list): Floats with the z-coordinates of the cell corners\n
        ncx (int): Number of cells in the x-dir\n
        ncz (int): Number of cells in the z-dir\n
        znr (list): Integers with the number of z-refinements per cell

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
