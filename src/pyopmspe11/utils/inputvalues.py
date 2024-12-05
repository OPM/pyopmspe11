# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""
Utiliy functions to set the requiried input values by pyopmspe11.
"""

import csv
import sys
from io import StringIO
from subprocess import PIPE, Popen
import numpy as np

try:
    import tomllib
except ImportError:
    pass


def process_input(dic, in_file):
    """
    Process the configuration file

    Args:
        dic (dict): Global dictionary\n
        in_file (str): Name of the input text file

    Returns:
        dic (dict): Modified global dictionary

    """
    if in_file.endswith(".toml"):
        if sys.version_info[1] < 11:
            print(
                "\nInput configuration files with toml extension requieres "
                + "a Python version of at least 3.11.\nYou could either use "
                + "configuration files with txt extension or use a higher "
                + "Python version.\nYour Python version is "
                + f"3.{sys.version_info[1]}.{sys.version_info[2]}.\n"
            )
            sys.exit()
        with open(in_file, "rb") as file:
            dic.update(tomllib.load(file))
        setcaseproperties(dic)
        postprocesstoml(dic)
    else:
        lol = []  # List of lines
        with open(in_file, "r", encoding="utf8") as file:
            for row in csv.reader(file, delimiter="#"):
                lol.append(row)
        readthefirstpart(lol, dic)
        setcaseproperties(dic)
        readthesecondpart(lol, dic)


def postprocesstoml(dic):
    """
    Perform the unit convertions and generation of needed variables

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    dic["noCells"] = [sum(dic["x_n"]), sum(dic["y_n"]), sum(dic["z_n"])]
    dic["diffusion"] = np.array(dic["diffusion"]) * 86400  # To [m^2/day]
    dic["noSands"] = len(dic["safu"])
    dic["tabdims"] = 1
    for row in dic["safu"]:
        dic["tabdims"] = max(dic["tabdims"], row[4])
    for i in range(len(dic["inj"])):  # To [s]
        dic["inj"][i][0] *= dic["time"]
        dic["inj"][i][1] *= dic["time"]
        dic["inj"][i][2] *= dic["time"]
    dic["wellCoord"][0][-1] = dic["dims"][2] - dic["wellCoord"][0][-1]
    dic["wellCoord"][1][-1] = dic["dims"][2] - dic["wellCoord"][1][-1]
    if dic["spe11"] == "spe11c":
        dic["wellCoordF"][0][-1] = dic["dims"][2] - dic["wellCoordF"][0][-1]
        dic["wellCoordF"][1][-1] = dic["dims"][2] - dic["wellCoordF"][1][-1]
    if dic["spe11"] in ["spe11b", "spe11c"]:
        dic["rockCond"] = np.array(dic["rockCond"]) * 86400.0 / 1e3  # To [kJ/(m day K)]


def setcaseproperties(dic):
    """
    Set the time scale, box, and sensor locations

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    dic["maxelevation"] = 0.0
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
        dic["maxelevation"] = 155.04166666666666  # at y = 2541 + 2/3
        corr = dic["elevation"] + 0.5 * dic["backElevation"]
        dic["sensors"] = [
            [4500.0, 2500.0, 655.0 - corr],
            [5100.0, 2500.0, 1255.0 - corr],
        ]
        dic["boxa"] = [[3300.0, 0.0, 0.0], [8300.0, 5000.0, 750.0]]
        dic["boxb"] = [[100.0, 0.0, 750.0], [3300.0, 5000.0, 1350.0]]
        dic["boxc"] = [[3300.0, 0.0, 250.0], [7800.0, 5000.0, 550.0]]
        dic["time"] = 31536000.0  # year to seconds


def readthefirstpart(lol, dic):
    """
    Process the lines from the flow executable and model parameters

    Args:
        lol (list): List of lines read from the input file\n
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    dic["flow"] = str(lol[1])[2:-2]  # Path to the flow executable
    row = (lol[4][0].strip()).split()
    dic["spe11"] = row[0]  # Name of the spe case (spe11a, spe11b, or spe11c)
    dic["version"] = row[1]  # OPM Flow version (release or master)
    row = (lol[5][0].strip()).split()
    dic["model"] = row[0]  # Model to run (immiscible, convective, or complete)
    dic["co2store"] = row[1]  # co2store implementation (gaswater or gasoil)
    dic["grid"] = str(lol[6][0]).strip()  # Grid (Cartesian, tensor, or corner-point)
    dic["dims"] = [float((lol[7][0].strip()).split()[j]) for j in range(3)]
    if dic["grid"] == "cartesian":
        dic["noCells"] = [int((lol[8 + j][0].strip()).split()[0]) for j in range(3)]
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
    dic["datum"] = float(row[0])  # Datum [m]
    dic["pressure"] = float(row[1])  # Pressure at the datum [Pa]
    dic["kzMult"] = float(row[2])  # Permeability multiplier in the z dir [-]
    row = list((lol[13][0].strip()).split())
    # Diffusion (in liquid and gas) [m^2/s] to [m^2/day]
    dic["diffusion"] = [float(row[0]) * 86400, float(row[1]) * 86400]
    if len(row) > 2:
        print(
            "Dispersion is now set per facie (Lines 37 to 43 in the configuration "
            + "file), not anymore in Line 14, e.g., see https://github.com/OPM/pyopm"
            + "spe11/blob/main/examples/finner_grids/spe11b.txt. Please update the "
            + "configuration file and try again to run pyopmspe11."
        )
        sys.exit()
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
    dic["index"] = 19  # Increase if more rows are added to the parameters part


def readthesecondpart(lol, dic):
    """
    Process the lines from the saturation functions until the end

    Args:
        lol (list): List of lines read from the input file\n
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

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
    for name in ["rock", "safu", "dispersion", "rockCond"]:
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
                str(float(row[1])),
                str(float(row[3])),
            ]
        )
        dic["dispersion"].append(float(row[5]))
        if dic["spe11"] != "spe11a":
            dic["rockCond"].append(float(row[7]) * 86400.0 / 1e3)
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
    dic["wellCoordF"] = columnf
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


def check_deck(dic):
    """
    Write unsupported features to the terminal if Flow from release is used

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    for value in dic["flow"].split():
        if "flow" in value:
            dic["only_flow"] = value
            break
    with Popen(
        args=f"{dic['only_flow']} --version", stdout=PIPE, shell=True
    ) as process:
        dic["flow_version"] = str(process.communicate()[0])[7:-3]
    if dic["flow_version"] == "2023.10":
        print(
            "\nYou are using Flow 2023.10. Please update to Flow 2024.04, or "
            + "build Flow from the master GitHub branches.\n"
        )
        sys.exit()
    if dic["flow_version"] == "2024.04" and dic["model"] == "convective":
        print(
            "\nThe convective option requires a Flow version newer than "
            + "22-08-2024. Please build Flow from the master GitHub branches.\n"
        )
        sys.exit()


def handle_tuning(dic):
    """
    If tuning is enable, then we write the entries corresponding to the maximum
    length of the next time step after event and the solver restart factor.

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    dic["tim_aft_eve"], dic["sol_res_fac"] = "", ""
    for value in dic["flow"].split():
        if "solver-restart-factor" in value:
            dic["sol_res_fac"] = float(value.split("=")[1])
            continue
        if "time-step-after-event-in-days" in value:
            dic["tim_aft_eve"] = float(value.split("=")[1])
