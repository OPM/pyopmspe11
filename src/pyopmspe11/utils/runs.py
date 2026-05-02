# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Utility functions for simulations, data processing, and plotting."""

import subprocess

from pyopmspe11.config.config import Config


def simulations(cfg: Config, flowfol: str) -> None:
    """Run OPM Flow."""
    data_file = f"{cfg.deckfol}/{cfg.fol.split('/')[-1].upper()}.DATA"
    flow_cmd = cfg.flow.split(" ") + [f"--output-dir={flowfol}", data_file]
    result = subprocess.run(flow_cmd, check=True)
    if result.returncode != 0:
        raise ValueError(f"Invalid result: {result.returncode}")


def plotting(cfg: Config) -> None:
    """Generate the figures."""
    plot_cmd = [
        "python3",
        f"{cfg.pat}/visualization/plotting.py",
        "-p",
        f"{cfg.fol}",
        "-d",
        f"{cfg.spe11}",
        "-g",
        f"{cfg.generate}",
        "-r",
        f"{cfg.resolution}",
        "-f",
        f"{cfg.subfolders}",
        "-t",
        f"{cfg.time_data}",
        "-n",
        "lower" if cfg.lower else "",
    ]
    print("\nPlot: Generation of png figures, please wait.")
    result = subprocess.run(plot_cmd, check=True)
    if result.returncode != 0:
        raise ValueError(f"Invalid result: {result.returncode}")


def data(cfg: Config) -> None:
    """Write the benchmark data."""
    data_cmd = [
        "python3",
        f"{cfg.pat}/visualization/data.py",
        "-p",
        f"{cfg.fol}",
        "-d",
        f"{cfg.spe11}",
        "-g",
        f"{cfg.generate}",
        "-r",
        f"{cfg.resolution}",
        "-t",
        f"{cfg.time_data}",
        "-w",
        f"{cfg.dt_data}",
        "-f",
        f"{cfg.subfolders}",
        "-n",
        "lower" if cfg.lower else "",
    ]
    print(
        "\nData: Generation of csv files following the SPE11 benchmark format, please wait."
    )
    result = subprocess.run(data_cmd, check=True)
    if result.returncode != 0:
        raise ValueError(f"Invalid result: {result.returncode}")
