# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""
Utiliy functions to set the requiried input values by pyopmcsp11.
"""

import csv
from io import StringIO
import numpy as np


def process_input(dic, in_file):
    """
    Function to process the input file

    Args:
        dic (dict): Global dictionary with required parameters
        in_file (str): Name of the input text file

    Returns:
        dic (dict): Global dictionary with new added parameters
    """
    lol = []  # List of lines
    with open(in_file, "r", encoding="utf8") as file:
        for row in csv.reader(file, delimiter="#"):
            lol.append(row)
    dic = readthefirstpart(lol, dic)
    dic = readthesecondpart(lol, dic)
    return dic


def readthefirstpart(lol, dic):
    """
    Function to process the lines from the flow executable and model parameters

    Args:
        lol (list): List of lines read from the input file
        dic (dict): Global dictionary with required parameters

    Returns:
        dic (dict): Global dictionary with new added parameters
    """
    dic["flow"] = str(lol[1])[2:-2]  # Path to the flow executable
    dic["csp11"] = (str(lol[4][0]).strip()).split()[
        0
    ]  # Name of the spe case (csp11a, csp11b, or csp11c)
    dic["version"] = (str(lol[4][0]).strip()).split()[
        1
    ]  # OPM Flow version (release or master)
    dic["model"] = str(lol[5][0]).strip()  # Model to run (immiscible or complete)
    dic["grid"] = str(lol[6][0]).strip()  # Type of grid (cartesian or corner-point)
    dic["dims"] = [float((lol[7][0].strip()).split()[j]) for j in range(3)]
    if dic["grid"] == "cartesian":
        dic["noCells"] = [int((lol[8 + j][0].strip()).split()[0]) for j in range(3)]
        dic["dsize"] = [1.0 * dic["dims"][i] / dic["noCells"][i] for i in range(3)]
    else:
        dic["x_n"] = np.genfromtxt(StringIO(lol[8][0]), delimiter=",", dtype=int)
        dic["y_n"] = np.genfromtxt(StringIO(lol[9][0]), delimiter=",", dtype=int)
        dic["z_n"] = np.genfromtxt(StringIO(lol[10][0]), delimiter=",", dtype=int)
        for ent in ["x_n", "y_n", "z_n"]:
            if dic[f"{ent}"].size == 1:
                dic[f"{ent}"] = [dic[f"{ent}"]]
        dic["noCells"] = [sum(dic["x_n"]), sum(dic["y_n"]), sum(dic["z_n"])]
    row = list((lol[11][0].strip()).split())
    # Temperature bottom and top rig [C]
    dic["temperature"] = [
        float(row[0]),
        float(row[1]),
    ]
    dic["pressure"] = float(lol[12][0])  # Pressure at the top [Pa]
    row = list((lol[13][0].strip()).split())
    # CO2 diffusion (in liquid and gas) [m^2/s]
    dic["diffusion"] = [float(row[0]), float(row[1])]
    row = list((lol[14][0].strip()).split())
    # Rock specific heat and density (for csp11b/c)
    dic["rockExtra"] = [float(row[0]), float(row[1])]
    row = list((lol[15][0].strip()).split())
    dic["pvAdded"] = float(row[0])  # Pore volume on lateral boundaries
    dic["widthBuffer"] = float(row[1])  # Width of the buffer cells [m]
    row = list((lol[16][0].strip()).split())
    dic["elevation"] = float(row[0])  # Elevation of the caprock [m]
    dic["backElevation"] = float(row[1])  # Elevation of back boundary [m]
    dic["noSands"] = 7  # No. saturation regions
    dic[
        "index"
    ] = 19  # Increase this if more rows are added to the model parameters part
    return dic


def readthesecondpart(lol, dic):
    """
    Function to process the lines from the saturation functions until the end

    Args:
        lol (list): List of lines read from the input file
        dic (dict): Global dictionary with required parameters

    Returns:
        dic (dict): Global dictionary with new added parameters
    """
    dic["krw"] = str(lol[dic["index"]][0])  # Wetting rel perm saturation function [-]
    dic["krn"] = str(
        lol[dic["index"] + 1][0]
    )  # Non-wetting rel perm saturation function [-]
    dic["pcap"] = str(
        lol[dic["index"] + 2][0]
    )  # Capillary pressure saturation function [Pa]
    dic["s_w"] = str(
        lol[dic["index"] + 3][0]
    )  # Capillary pressure saturation function [Pa]
    dic["index"] += 7
    for name in ["rock", "safu", "rockCond"]:
        dic[name] = []
    dic["tabdims"] = 1
    for i in range(dic["noSands"]):  # Saturation function values
        row = list((lol[dic["index"] + i][0].strip()).split())
        dic["safu"].append(
            [
                float(row[1]),
                float(row[3]),
                float(row[5]),
                float(row[7]),
                int(row[9]),
            ]
        )
        dic["tabdims"] = max(dic["tabdims"], int(row[9]))
    dic["index"] += 3 + dic["noSands"]
    for i in range(dic["noSands"]):  # Rock values
        row = list((lol[dic["index"] + i][0].strip()).split())
        dic["rock"].append(
            [
                row[1],
                row[3],
            ]
        )
        if dic["csp11"] != "csp11a":
            dic["rockCond"].append(
                [
                    float(row[5]) * 86400.0 / 1e3,
                ]
            )

    dic["index"] += 3 + dic["noSands"]
    column = []
    columnf = []
    dic["radius"] = []
    for i in range(len(lol) - dic["index"]):
        if not lol[dic["index"] + i]:
            break
        row = list((lol[dic["index"] + i][0].strip()).split())
        dic["radius"].append(float(row[0]))
        column.append(
            [
                float(row[1]),
                float(row[2]),
                dic["dims"][2] - float(row[3]),
            ]
        )
        if dic["csp11"] == "csp11c":
            columnf.append(
                [
                    float(row[4]),
                    float(row[5]),
                    dic["dims"][2] - float(row[6]),
                ]
            )
    dic["wellCoord"] = column
    dic["wellCoordf"] = columnf
    dic["index"] += len(dic["wellCoord"]) + 3
    column = []
    for i in range(len(lol) - dic["index"]):
        if not lol[dic["index"] + i]:
            break
        row = list((lol[dic["index"] + i][0].strip()).split())
        column.append([float(row[j]) for j in range(3 + 3 * len(dic["wellCoord"]))])
    dic["inj"] = column
    return dic
