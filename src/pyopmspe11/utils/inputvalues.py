# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT
# pylint: disable=R0915

"""
Utility functions to set required input values for pyopmspe11.
"""

import csv
import sys
import os
from io import StringIO
from typing import Any
import subprocess
from subprocess import PIPE
import numpy as np

from pyopmspe11.config.config import Config

try:
    import tomllib
except ImportError:
    pass


def process_input(cli: dict) -> Config:
    """Process configuration input."""
    msg1 = (
        "\nAfter the pyopmspe11 2025.04 release, the CO2STORE functionality only "
        + "uses the gaswater implementation, not the gasoil implementation.\nThen "
        + "either remove the gasoil text in your configuration file, or use an "
        + "older release of pyopmspe11.\nThe execution of pyopmspe11 will continue "
        + "setting the deck with the gaswater implementation."
    )
    msg2 = (
        "\nAfter the pyopmspe11 2025.04 release, column 3 for the maximum solver time "
        + "step in the injection has been moved to the end of the column, including the "
        + "items for the TUNING keyword, which gives more control when setting "
        + "the simulations. Please see the configuration files in the examples and "
        + "online documentation, and update your configuration file accordingly.\n"
    )
    if cli["input"].endswith(".toml"):
        with open(cli["input"], "rb") as f:
            cfg_file = tomllib.load(f)
    else:
        lines = []
        with open(cli["input"], "r", encoding="utf8") as file:
            for row in csv.reader(file, delimiter="#"):
                lines.append(row)
        cfg_file = load_config_txt(lines)
    cfg = Config(
        fol=os.path.abspath(cli["output"]),
        generate=cli["generate"],
        mode=cli["mode"],
        resolution=cli["resolution"],
        time_data=cli["time"],
        dt_data=float(cli["write"]),
        lower=cli["neighbourhood"],
        subfolders=cli["subfolders"],
        **cfg_file,
    )
    time = setcaseproperties(cfg)
    postprocesstoml(cfg, time, msg1, msg2)
    for value in cfg.flow.split():
        if "--enable-tuning" in value and value[16:] in ["true", "True", "1"]:
            cfg.tuning = True
            break

    return cfg


def load_config_txt(lines: list) -> dict:
    """Read the txt  (the types are checked when creating the Config data class)"""
    dic: dict[str, Any] = {"flow": str(lines[1])[2:-2]}
    row = lines[4][0].strip().split()
    dic["spe11"] = row[0]
    dic["version"] = row[1]
    row = lines[5][0].strip().split()
    dic["model"] = row[0]
    dic["grid"] = lines[6][0].strip()
    split7 = lines[7][0].strip().split()
    dic["dims"] = [float(split7[0]), float(split7[1]), float(split7[2])]
    dic["x_n"] = np.genfromtxt(StringIO(lines[8][0]), delimiter=",", dtype=int)
    dic["y_n"] = np.genfromtxt(StringIO(lines[9][0]), delimiter=",", dtype=int)
    dic["z_n"] = np.genfromtxt(StringIO(lines[10][0]), delimiter=",", dtype=int)
    for key in ["x_n", "y_n", "z_n"]:
        if np.size(dic[key]) == 1:
            dic[key] = [int(dic[key])]
    row = lines[11][0].strip().split()
    dic["temperature"] = [float(row[0]), float(row[1])]
    row = lines[12][0].strip().split()
    dic["datum"] = float(row[0])
    dic["pressure"] = float(row[1])
    dic["kzMult"] = float(row[2])
    row = lines[13][0].strip().split()
    dic["diffusion"] = [float(row[0]), float(row[1])]
    row = lines[14][0].strip().split()
    dic["rockExtra"] = [float(row[0]), float(row[1])]
    row = lines[15][0].strip().split()
    dic["spe11aBC"] = float(row[0])
    dic["pvAdded"] = float(row[1])
    dic["widthBuffer"] = float(row[2])
    row = lines[16][0].strip().split()
    dic["elevation"] = float(row[0])
    dic["backElevation"] = float(row[1])
    idx = 19
    dic["krw"] = str(lines[idx][0])
    dic["krn"] = str(lines[idx + 1][0])
    dic["pcap"] = str(lines[idx + 2][0])
    dic["s_w"] = str(lines[idx + 3][0])
    idx += 7
    dic["rock"] = [[] for _ in range(7)]
    dic["safu"] = [[] for _ in range(7)]
    dic["dispersion"] = [0.0 for _ in range(7)]
    dic["rockCond"] = [0.0 for _ in range(7)]
    for ind in range(7):
        row = lines[idx + ind][0].strip().split()
        dic["safu"][ind] = [
            float(row[1]),
            float(row[3]),
            float(row[5]),
            float(row[7]),
            int(row[9]),
        ]
    idx += 10
    for ind in range(7):
        row = lines[idx + ind][0].strip().split()
        dic["rock"][ind] = [float(row[1]), float(row[3])]
        dic["dispersion"][ind] = float(row[5])
        if dic["spe11"] != "spe11a":
            dic["rockCond"][ind] = float(row[7])
    idx += 10
    dic["wellCoord"], dic["wellCoordF"], dic["radius"] = [], [], []
    for ind in range(len(lines) - idx):
        if not lines[idx + ind]:
            break
        row = lines[idx + ind][0].strip().split()
        dic["radius"].append(float(row[0]))
        dic["wellCoord"].append([float(row[1]), float(row[2]), float(row[3])])
        if dic["spe11"] == "spe11c":
            dic["wellCoordF"].append([float(row[4]), float(row[5]), float(row[6])])
    idx += len(dic["wellCoord"]) + 3
    injections, tunning = [], []
    for ind in range(len(lines) - idx):
        if not lines[idx + ind]:
            break
        row = lines[idx + ind][0].strip().split()
        entry = [float(row[0]), float(row[1])] + [float(row[j]) for j in range(2, 8)]
        if len(row) > 8:
            parts = " ".join(row[8:]).split("/")
            for val in parts:
                tunning.append(val.strip().replace("'", "").replace('"', ""))
        injections.append(entry + tunning)
    dic["inj"] = injections
    return dic


def postprocesstoml(cfg: Config, time: float, msg1: str, msg2: str) -> None:
    """Convert units and generate variables."""
    cfg.nxyz = [sum(cfg.x_n), sum(cfg.y_n), sum(cfg.z_n)]
    cfg.diffusion = [val * 86400 for val in cfg.diffusion]
    for inj in cfg.inj:
        inj[0] *= time
        inj[1] *= time
    dims_z = cfg.dims[2]
    cfg.wellCoord[0][-1] = dims_z - cfg.wellCoord[0][-1]
    cfg.wellCoord[1][-1] = dims_z - cfg.wellCoord[1][-1]
    if cfg.spe11 == "spe11c":
        assert cfg.wellCoordF is not None
        cfg.wellCoordF[0][-1] = dims_z - cfg.wellCoordF[0][-1]
        cfg.wellCoordF[1][-1] = dims_z - cfg.wellCoordF[1][-1]
    if not hasattr(cfg, "rockCond"):
        cfg.rockCond = [0.0] * 7
    if not hasattr(cfg, "rockExtra"):
        cfg.rockExtra = [0.0, 0.0]
    if cfg.rockCond:
        cfg.rockCond = [val * 86400.0 / 1e3 for val in cfg.rockCond]
    if getattr(cfg, "co2store", None) == "gasoil":
        print(msg1)
    for inj in cfg.inj:
        if len(inj) == 9 and not isinstance(inj[-1], str):
            print(msg2)
            sys.exit()
        if len(inj) >= 9 and isinstance(inj[-1], str):
            parts = inj[-1].split("/")
            inj[-1] = parts[0].strip()
            for extra in parts[1:]:
                inj.append(extra.strip())


def setcaseproperties(cfg: Config) -> float:
    """Set time scale, boxes, and sensors."""
    mult = 0.995 if cfg.lower and cfg.grid != "corner-point" else 1
    if cfg.spe11 == "spe11a":
        cfg.sensors = [[1.5, 0.005, mult * 0.5], [1.7, 0.005, 1.1]]
        cfg.boxa = [[1.1, 0.0, 0.0], [2.8, 0.01, 0.6]]
        cfg.boxb = [[0.0, 0.0, 0.6], [1.1, 0.01, 1.2]]
        cfg.boxc = [[1.1, 0.0, 0.1], [2.6, 0.01, 0.4]]
        time = 3600.0
    elif cfg.spe11 == "spe11b":
        cfg.sensors = [[4500.0, 0.5, mult * 500], [5100.0, 0.5, 1100.0]]
        cfg.boxa = [[3300.0, 0.0, 0.0], [8300.0, 1.0, 600.0]]
        cfg.boxb = [[100.0, 0.0, 600.0], [3300.0, 1.0, 1200.0]]
        cfg.boxc = [[3300.0, 0.0, 100.0], [7800.0, 1.0, 400.0]]
        time = 31536000.0
    else:
        cfg.maxelevation = 155.04166666666666
        assert cfg.elevation is not None
        assert cfg.backElevation is not None
        correction = cfg.elevation + 0.5 * cfg.backElevation
        cfg.sensors = [
            [4500.0, 2500.0, mult * (655.0 - correction)],
            [5100.0, 2500.0, 1255.0 - correction],
        ]
        cfg.boxa = [[3300.0, 0.0, 0.0], [8300.0, 5000.0, 750.0]]
        cfg.boxb = [[100.0, 0.0, 750.0], [3300.0, 5000.0, 1350.0]]
        cfg.boxc = [[3300.0, 0.0, 250.0], [7800.0, 5000.0, 550.0]]
        time = 31536000.0
    return time


def check_deck(cfg: Config) -> None:
    """Check Flow release compatibility."""
    only_flow = ""
    for value in cfg.flow.split():
        if "flow" in value:
            only_flow = value
            break
    result = subprocess.run(
        f"{only_flow} --version", stdout=PIPE, shell=True, check=False
    )
    flow_version = result.stdout.decode(errors="ignore")[7:-3]
    for forbidden in ["2025.04", "2024.10", "2024.04"]:
        if flow_version == forbidden:
            print(
                f"\nYou are using Flow {forbidden}. Please update to Flow 2025.10, "
                + "or build Flow from the master GitHub branches.\n"
            )
            sys.exit()
