# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""
Utiliy functions to set the requiried input values by pyopmspe11.
"""

import csv
from io import StringIO
from subprocess import PIPE, Popen
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
    if dic["spe11"] == "spe11a":
        dic["sensors"] = [[1.5, 0.005, 0.5], [1.7, 0.005, 1.1]]
        dic["boxa"] = [[1.1, 0.0, 0.0], [2.8, 0.01, 0.6]]
        dic["boxb"] = [[0.0, 0.0, 0.6], [1.1, 0.01, 1.2]]
        dic["boxc"] = [[1.1, 0.0, 0.1], [2.6, 0.01, 0.4]]
        dic["time"] = 3600.0  # hour to seconds
    elif dic["spe11"] == "spe11b":
        dic["sensors"] = [[4500.0, 0.5, 500], [5100.0, 0.5, 1100.0]]
        dic["boxa"] = [[3300.0, 0.0, 0.0], [8300.0, 1.0, 600.0]]
        dic["boxb"] = [[100.0, 0.0, 600.0], [3300.0, 1.0, 1200.0]]
        dic["boxc"] = [[3300.0, 0.0, 100.0], [7800.0, 1.0, 400.0]]
        dic["time"] = 31536000.0  # year to seconds
    else:
        corr = dic["elevation"] + 0.5 * dic["backElevation"]
        dic["sensors"] = [
            [4500.0, 2500.0, 655.0 - corr],
            [5100.0, 2500.0, 1255.0 - corr],
        ]
        dic["boxa"] = [[3300.0, 0.0, 0.0], [8300.0, 5000.0, 750.0]]
        dic["boxb"] = [[100.0, 0.0, 750.0], [3300.0, 5000.0, 1350.0]]
        dic["boxc"] = [[3300.0, 0.0, 250.0], [7800.0, 5000.0, 550.0]]
        dic["time"] = 31536000.0  # year to seconds
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
    row = (lol[4][0].strip()).split()
    dic["spe11"] = row[0]  # Name of the spe case (spe11a, spe11b, or spe11c)
    dic["version"] = row[1]  # OPM Flow version (release or master)
    row = (lol[5][0].strip()).split()
    dic["model"] = row[0]  # Model to run (immiscible or complete)
    dic["co2store"] = row[1]  # co2store implementation (gaswater or gasoil)
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
    row = list((lol[12][0].strip()).split())
    dic["pressure"] = float(row[0])  # Pressure at the top [Pa]
    dic["kzMult"] = float(row[1])  # Permeability multiplier in the z dir [-]
    row = list((lol[13][0].strip()).split())
    # Diffusion (in liquid and gas) [m^2/s] to [m^2/day]
    dic["diffusion"] = [float(row[0]) * 86400, float(row[1]) * 86400]
    dic["dispersion"] = float(row[2])  # Dispersion [m]
    row = list((lol[14][0].strip()).split())
    # Rock specific heat and density (for spe11b/c)
    dic["rockExtra"] = [float(row[0]), float(row[1])]
    row = list((lol[15][0].strip()).split())
    dic["spe11aBC"] = float(row[0])  # Boundary on top (spe11a [free or added pv])
    dic["pvAdded"] = float(row[1])  # Pore volume on lateral boundaries
    dic["widthBuffer"] = float(row[2])  # Width of the buffer cells [m]
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
        if dic["spe11"] != "spe11a":
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
        if dic["spe11"] == "spe11c":
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
        column.append(
            [
                float(row[0]) * dic["time"],
                float(row[1]) * dic["time"],
                float(row[2]) * dic["time"],
            ]
            + [float(row[j]) for j in range(3, 3 + 3 * len(dic["wellCoord"]))]
        )
    dic["inj"] = column
    return dic


def check_deck(dic):
    """Write unsupported features to the terminal if flow from release is used"""
    for value in dic["flow"].split():
        if "flow" in value:
            flow = value
            break
    with Popen(args=f"{flow} --version", stdout=PIPE, shell=True) as process:
        dic["flow_version"] = str(process.communicate()[0])[7:-3]
    if dic["flow_version"] == "2023.10":
        if dic["co2store"] == "gaswater" and dic["spe11"] != "spe11a":
            print(
                "\nDiffusion is not supported for gaswater + energy systems in "
                + "flow 2023.10.\nThen diffusion is not included in the generated "
                + "deck.\nEither select the co2store gasoil implementation or build "
                + "flow from the master GitHub branches.\n"
            )
        if dic["dispersion"] > 0:
            print(
                "\nDispersion is not supported in flow 2023.10.\nThen dispersion is "
                + "not included in the generated deck.\nBuild flow from the master "
                + "GitHub branches to include dispersion in the simulations.\n"
            )
