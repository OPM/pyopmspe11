# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""
Utiliy function for the grid and finding the wells i,j, and k ids.
"""

import csv
import numpy as np
import pandas as pd
from resdata.grid import Grid
from shapely.geometry import Polygon


def grid(dic):
    """
    Function to handle grid types

    Args:
        dic (dict): Global dictionary with required parameters

    Returns:
        dic (dict): Global dictionary with new added parameters

    """
    if dic["grid"] == "corner-point":
        dic = corner(dic)
    elif dic["grid"] == "cartesian":
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
        for name, size in zip(["xmx", "ymy", "zmz"], ["dx", "dy", "dz"]):
            dic[f"{name}_center"] = (dic[f"{name}"][1:] + dic[f"{name}"][:-1]) / 2.0
            dic[f"{size}"] = dic[f"{name}"][1:] - dic[f"{name}"][:-1]

    return dic


def structured_handling_spe11a(dic):
    """
    Function to locate sand and well positions in the tensor/cartesian grid

    Args:
        dic (dict): Global dictionary with required parameters

    Returns:
        dic (dict): Global dictionary with new added parameters

    """
    sensor1, sensor2 = [], []
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
            dic = boxes(dic, dic["xmx_center"][i], dic["zmz_center"][k], i)
            dic["satnum"].append(dic["ids_gmsh"][fgl][0])
            dic["permx"].append(dic["rock"][int(dic["ids_gmsh"][fgl][0]) - 1][0])
            dic["poro"].append(dic["rock"][int(dic["ids_gmsh"][fgl][0]) - 1][1])
    dic["fipnum"][pd.Series(sensor1).argmin()] = "5"
    dic["fipnum"][pd.Series(sensor2).argmin()] = "6"
    dic = wells(dic)
    return dic


def structured_handling_spe11bc(dic):
    """
    Function to locate sand and well positions in the tensor/cartesian grid

    Args:
        dic (dict): Global dictionary with required parameters

    Returns:
        dic (dict): Global dictionary with new added parameters

    """
    sensor1, sensor2 = [], []
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
                z_c += map_z(dic, 0)
            dic = boxes(dic, dic["xmx_center"][i], z_c, i)
            dic["satnum"].append(dic["ids_gmsh"][fgl][0])
            dic["permx"].append(dic["rock"][int(dic["ids_gmsh"][fgl][0]) - 1][0])
            dic["thconr"].append(
                f"{dic['rockCond'][int(dic['ids_gmsh'][fgl][0])-1][0]}"
            )
            if i == 0 and (
                int(dic["ids_gmsh"][fgl][0]) != 1 and int(dic["ids_gmsh"][fgl][0]) != 7
            ):
                dic["poro"].append(f"{dic['pvAdded']*dic['ymy'][1]*dic['dz'][k] :E}")
            elif i == dic["noCells"][0] - 1 and (
                int(dic["ids_gmsh"][fgl][0]) != 1 and int(dic["ids_gmsh"][fgl][0]) != 7
            ):
                dic["poro"].append(f"{dic['pvAdded']*dic['ymy'][1]*dic['dz'][k] :E}")
            else:
                dic["poro"].append(dic["rock"][int(dic["ids_gmsh"][fgl][0]) - 1][1])
        for j in range(dic["noCells"][1] - 1):
            for names in ["satnum", "poro", "permx", "thconr"]:
                dic[f"{names}"].extend(dic[f"{names}"][-dic["noCells"][0] :])
            for i_i in range(dic["noCells"][0]):
                sensor1.append(
                    (dic["xmx_center"][i_i] - dic["sensors"][0][0]) ** 2
                    + (dic["ymy_center"][j] - dic["sensors"][0][1]) ** 2
                    + (dic["zmz_center"][k] + dic["sensors"][0][2] - dic["dims"][2])
                    ** 2
                )
                sensor2.append(
                    (dic["xmx_center"][i_i] - dic["sensors"][1][0]) ** 2
                    + (dic["ymy_center"][j] - dic["sensors"][1][1]) ** 2
                    + (dic["zmz_center"][k] + dic["sensors"][1][2] - dic["dims"][2])
                    ** 2
                )
                z_c = dic["zmz_center"][k]
                if dic["spe11"] == "spe11c":
                    z_c += map_z(dic, j + 1)
                dic = boxes(dic, dic["xmx_center"][i_i], z_c, i_i)
    dic["fipnum"][pd.Series(sensor1).argmin()] = "5"
    dic["fipnum"][pd.Series(sensor2).argmin()] = "6"
    dic = wells(dic)
    return dic


def corner_point_handling_spe11a(dic):
    """
    Function to locate sand and well positions in a corner-point grid

    Args:
        dic (dict): Global dictionary with required parameters

    Returns:
        dic (dict): Global dictionary with new added parameters

    """
    well1, well2, sensor1, sensor2 = [], [], [], []
    gridf = Grid(f"{dic['exe']}/{dic['fol']}/flow/INITIAL.EGRID")
    dic["wellijk"] = [[] for _ in range(len(dic["wellCoord"]))]
    for cell in gridf.cells():
        fgl = pd.Series(
            (
                (dic["cxc1"] - cell.coordinate[0]) ** 2
                + (dic["czc1"] - cell.coordinate[2]) ** 2
            )
        ).argmin()
        well1.append(
            (dic["wellCoord"][0][0] - cell.coordinate[0]) ** 2
            + (dic["wellCoord"][0][2] - cell.coordinate[2]) ** 2
        )
        well2.append(
            (dic["wellCoord"][1][0] - cell.coordinate[0]) ** 2
            + (dic["wellCoord"][1][2] - cell.coordinate[2]) ** 2
        )
        sensor1.append(
            (cell.coordinate[0] - dic["sensors"][0][0]) ** 2
            + (cell.coordinate[2] + dic["sensors"][0][2] - dic["dims"][2]) ** 2
        )
        sensor2.append(
            (cell.coordinate[0] - dic["sensors"][1][0]) ** 2
            + (cell.coordinate[2] + dic["sensors"][1][2] - dic["dims"][2]) ** 2
        )
        dic = boxes(dic, cell.coordinate[0], cell.coordinate[2], cell.i)
        dic["satnum"].append(dic["ids_gmsh"][fgl][0])
        dic["permx"].append(dic["rock"][int(dic["ids_gmsh"][fgl][0]) - 1][0])
        dic["poro"].append(dic["rock"][int(dic["ids_gmsh"][fgl][0]) - 1][1])
    idwell1 = pd.Series(well1).argmin()
    idwell2 = pd.Series(well2).argmin()
    dic["wellijk"][0] = [gridf.cell(idwell1).i + 1, 1, gridf.cell(idwell1).k + 1]
    dic["wellijk"][1] = [gridf.cell(idwell2).i + 1, 1, gridf.cell(idwell2).k + 1]
    return dic


def corner_point_handling_spe11bc(dic):
    """
    Function to locate sand and well positions in a corner-point grid

    Args:
        dic (dict): Global dictionary with required parameters

    Returns:
        dic (dict): Global dictionary with new added parameters

    """
    well1, well2, sensor1, sensor2, xtemp = [], [], [], [], []
    gridf = Grid(f"{dic['exe']}/{dic['fol']}/flow/INITIAL.EGRID")
    dic["wellijk"] = [[] for _ in range(len(dic["wellCoord"]))]
    for cell in gridf.cells():
        xtemp.append(cell.coordinate[0])
        fgl = pd.Series(
            (
                (dic["cxc1"] - cell.coordinate[0]) ** 2
                + (dic["czc1"] - cell.coordinate[2]) ** 2
            )
        ).argmin()
        well1.append(
            (dic["wellCoord"][0][0] - cell.coordinate[0]) ** 2
            + (dic["wellCoord"][0][2] - cell.coordinate[2]) ** 2
        )
        well2.append(
            (dic["wellCoord"][1][0] - cell.coordinate[0]) ** 2
            + (dic["wellCoord"][1][2] - cell.coordinate[2]) ** 2
        )
        sensor1.append(
            (cell.coordinate[0] - dic["sensors"][0][0]) ** 2
            + (dic["ymy_center"][0] - dic["sensors"][0][1]) ** 2
            + (cell.coordinate[2] + dic["sensors"][0][2] - dic["dims"][2]) ** 2
        )
        sensor2.append(
            (cell.coordinate[0] - dic["sensors"][1][0]) ** 2
            + (dic["ymy_center"][0] - dic["sensors"][1][1]) ** 2
            + (cell.coordinate[2] + dic["sensors"][1][2] - dic["dims"][2]) ** 2
        )
        z_c = cell.coordinate[2]
        if dic["spe11"] == "spe11c":
            z_c += map_z(dic, cell.j)
        dic = boxes(dic, cell.coordinate[0], z_c, cell.i)
        dic["satnum"].append(dic["ids_gmsh"][fgl][0])
        dic["permx"].append(dic["rock"][int(dic["ids_gmsh"][fgl][0]) - 1][0])
        dic["thconr"].append(f"{dic['rockCond'][int(dic['ids_gmsh'][fgl][0])-1][0]}")
        if cell.i == 0 and (
            int(dic["ids_gmsh"][fgl][0]) != 1 and int(dic["ids_gmsh"][fgl][0]) != 7
        ):
            dic["poro"].append(
                f"{dic['pvAdded']*cell.dimension[1]*cell.dimension[2] :E}"
            )
        elif cell.i == dic["noCells"][0] - 1 and (
            int(dic["ids_gmsh"][fgl][0]) != 1 and int(dic["ids_gmsh"][fgl][0]) != 7
        ):
            dic["poro"].append(
                f"{dic['pvAdded']*cell.dimension[1]*cell.dimension[2] :E}"
            )
        else:
            dic["poro"].append(dic["rock"][int(dic["ids_gmsh"][fgl][0]) - 1][1])
        if cell.i > 0 and cell.i == dic["noCells"][0] - 1:
            for j in range(dic["noCells"][1] - 1):
                for names in ["satnum", "poro", "permx", "thconr"]:
                    dic[f"{names}"].extend(dic[f"{names}"][-dic["noCells"][0] :])
                for i_i in range(dic["noCells"][0]):
                    sensor1.append(
                        (xtemp[i_i] - dic["sensors"][0][0]) ** 2
                        + (dic["ymy_center"][j + 1] - dic["sensors"][0][1]) ** 2
                        + (cell.coordinate[2] + dic["sensors"][0][2] - dic["dims"][2])
                        ** 2
                    )
                    sensor2.append(
                        (xtemp[i_i] - dic["sensors"][1][0]) ** 2
                        + (dic["ymy_center"][j + 1] - dic["sensors"][1][1]) ** 2
                        + (cell.coordinate[2] + dic["sensors"][1][2] - dic["dims"][2])
                        ** 2
                    )
                    z_c = cell.coordinate[2]
                    if dic["spe11"] == "spe11c":
                        z_c += map_z(dic, j + 1)
                    dic = boxes(dic, xtemp[i_i], z_c, i_i)
            xtemp = []
    well1 = pd.Series(well1).argmin()
    well2 = pd.Series(well2).argmin()
    dic["fipnum"][pd.Series(sensor1).argmin()] = "5"
    dic["fipnum"][pd.Series(sensor2).argmin()] = "6"
    dic["wellijk"][0] = [gridf.cell(well1).i + 1, 1, gridf.cell(well1).k + 1]
    dic["wellijk"][1] = [gridf.cell(well2).i + 1, 1, gridf.cell(well2).k + 1]
    # Work in process to implement properly this for the corner-point grid in spe11c
    if dic["spe11"] == "spe11c":
        dic["wellijkf"] = [[] for _ in range(len(dic["wellCoord"]))]
        dic["wellijkf"][0] = [dic["wellijk"][0][0], 1, dic["wellijk"][0][2]]
        dic["wellijkf"][1] = [dic["wellijk"][1][0], 1, dic["wellijk"][1][2]]
        dic["ymy_center"] = (np.array(dic["ymy"][1:]) + np.array(dic["ymy"][:-1])) / 2.0
        wji = pd.Series(abs(dic["wellCoord"][0][1] - dic["ymy_center"])).argmin() + 1
        wjf = pd.Series(abs(dic["wellCoordf"][0][1] - dic["ymy_center"])).argmin() + 1
        dic["wellijk"][0][1] = wji
        dic["wellijk"][1][1] = wji
        dic["wellijkf"][0][1] = wjf
        dic["wellijkf"][1][1] = wjf
        dic["wellkh"] = [
            dic["wellijk"][1][2]
            for _ in range(dic["wellijk"][1][1], dic["wellijkf"][1][1] + 1)
        ]
    return dic


def boxes(dic, x_c, z_c, idx):
    """Find the box positions"""
    if (
        (dic["dims"][2] - z_c >= dic["boxb"][0][2])
        & (dic["dims"][2] - z_c <= dic["boxb"][1][2])
        & (x_c >= dic["boxb"][0][0])
        & (x_c <= dic["boxb"][1][0])
    ):
        dic["fipnum"].append("3")
    elif (
        (dic["dims"][2] - z_c >= dic["boxc"][0][2])
        & (dic["dims"][2] - z_c <= dic["boxc"][1][2])
        & (x_c >= dic["boxc"][0][0])
        & (x_c <= dic["boxc"][1][0])
    ):
        dic["fipnum"].append("4")
    elif (
        (dic["dims"][2] - z_c >= dic["boxa"][0][2])
        & (dic["dims"][2] - z_c <= dic["boxa"][1][2])
        & (x_c >= dic["boxa"][0][0])
        & (x_c <= dic["boxa"][1][0])
    ):
        dic["fipnum"].append("2")
    elif dic["spe11"] != "spe11a" and idx in (0, dic["noCells"][0] - 1):
        dic["fipnum"].append("7")
    else:
        dic["fipnum"].append("1")
    return dic


def positions(dic):
    """
    Function to locate sand and well positions

    Args:
        dic (dict): Global dictionary with required parameters

    Returns:
        dic (dict): Global dictionary with new added parameters

    """
    dic = getfacies(dic)
    for names in ["satnum", "poro", "permx", "thconr", "fipnum"]:
        dic[f"{names}"] = []
    if dic["grid"] == "corner-point":
        if dic["spe11"] == "spe11a":
            dic = corner_point_handling_spe11a(dic)
        else:
            dic = corner_point_handling_spe11bc(dic)
    else:
        if dic["spe11"] == "spe11a":
            dic = structured_handling_spe11a(dic)
        else:
            dic = structured_handling_spe11bc(dic)
    return dic


def wells(dic):
    """
    Function to find the wells index

    Args:
        dic (dict): Global dictionary with required parameters

    Returns:
        dic (dict): Global dictionary with new added parameters

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
                        abs(dic["wellCoordf"][j][k] - dic[f"{axis}_center"])
                    ).argmin()
                    + 1
                )
            well_coord = dic["wellCoord"][j][2]
            midpoints = dic["zmz_center"] + map_z(dic, dic["wellijk"][j][1])
            dic["wellijk"][j].append(
                pd.Series(
                    abs(
                        well_coord
                        + dic["elevation"]
                        + 0.5 * dic["backElevation"]
                        - midpoints
                    )
                ).argmin()
                + 1
            )
    if dic["spe11"] == "spe11c":
        dic["wellkh"] = []
        for j in range(dic["wellijk"][1][1], dic["wellijkf"][1][1] + 1):
            midpoints = dic["zmz_center"] + map_z(dic, j)
            dic["wellkh"].append(
                pd.Series(
                    abs(
                        dic["wellCoord"][1][2]
                        + map_z(dic, dic["wellijk"][1][1])
                        - midpoints
                    )
                ).argmin()
                + 1
            )
    return dic


def map_z(dic, j):
    """
    Function to return the z position of the parabola

    Args:
        dic (dict): Global dictionary with required parameters
        j : cell id along the y axis

    Returns:
        z: Position

    """
    z_pos = (
        dic["elevation"]
        - dic["elevation"]
        * (1.0 - (dic["ymy_center"][j - 1] / (0.5 * dic["dims"][1]) - 1) ** 2.0)
        - dic["ymy_center"][j - 1] * dic["backElevation"] / dic["dims"][1]
    )
    return z_pos


def getfacies(dic):
    """
    Function to read the reference gmsh file

    Args:
        dic (dict): Global dictionary with required parameters

    Returns:
        dic (dict): Global dictionary with new added parameters

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
    return dic


def get_lines(dic):
    """
    Function to read the points in the z-surface lines

    Args:
        dic (dict): Global dictionary with required parameters

    Returns:
        lines: List with the coordinates of the lines

    """
    with open(
        f"{dic['pat']}/reference_mesh/lines_coordinates.geo",
        "r",
        encoding="utf8",
    ) as file:
        lol = file.readlines()
    lines = []
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
    Function to create a spe11 corner-point grid

    Args:
        dic (dict): Global dictionary with required parameters

    Returns:
        dic (dict): Global dictionary with new added parameters

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
    dic["noCells"][1] = len(dic["ymy"]) - 1
    if (dic["spe11"] == "spe11b" or dic["spe11"] == "spe11c") and 1.1 * dic[
        "widthBuffer"
    ] < dic["xmx"][1]:
        dic["xmx"] = np.insert(dic["xmx"], 1, dic["widthBuffer"])
        dic["xmx"] = np.insert(
            dic["xmx"], len(dic["xmx"]) - 1, dic["xmx"][-1] - dic["widthBuffer"]
        )
    for xcor in dic["xmx"]:
        for i, lcor in enumerate(lines):
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
    return dic


def refinement_z(xci, zci, ncx, ncz, znr):
    """
    Refinment of the grid in the z-dir

    Args:
        xci (list): List of floats with the x-coordinates of the cell corners
        zci (list): List of floats with the z-coordinates of the cell corners
        ncx (int): number of cells in the x-dir
        ncz (int): number of cells in the z-dir
        znr (list): List of integers with the number of z-refinments per cell

    Returns:
        xcr (list): List of floats with the new x-coordinates of the cell corners
        zcr (list): List of floats with the new z-coordinates of the cell corners
        ncx (int): new number of cells in the x-dir
        ncz (int): new number of cells in the z-dir
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
