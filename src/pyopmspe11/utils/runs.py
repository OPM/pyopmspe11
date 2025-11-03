# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""
Utiliy functions for the simulations, data processing, and plotting.
"""
import os
import subprocess


def simulations(dic):
    """
    Run OPM Flow

    Args:
        dic (dict): Global dictionary

    Returns:
        None

    """
    os.system(
        f"{dic['flow']} --output-dir={dic['flowf']} "
        f"{dic['deckf']}/{dic['fol'].split('/')[-1].upper()}.DATA & wait\n"
    )


def plotting(dic):
    """
    Generate the figures

    Args:
        dic (dict): Global dictionary

    Returns:
        None

    """
    plot_exe = [
        "python3",
        f"{dic['pat']}/visualization/plotting.py",
        "-p " + f"{dic['fol']}",
        "-d " + f"{dic['spe11']}",
        "-g " + f"{dic['generate']}",
        "-r " + f"{dic['resolution']}",
        "-s " + f"{dic['showpywarn']}",
        "-l " + f"{dic['latex']}",
        "-f " + f"{dic['subfolders']}",
        "-t " + f"{dic['time_data']}",
    ]
    print("\nPlot: Generation of png figures, please wait.")
    prosc = subprocess.run(plot_exe, check=True)
    if prosc.returncode != 0:
        raise ValueError(f"Invalid result: { prosc.returncode }")


def data(dic):
    """
    Write the benchmark data

    Args:
        dic (dict): Global dictionary

    Returns:
        None

    """
    data_exe = [
        "python3",
        f"{dic['pat']}/visualization/data.py",
        "-p " + f"{dic['fol']}",
        "-d " + f"{dic['spe11']}",
        "-g " + f"{dic['generate']}",
        "-r " + f"{dic['resolution']}",
        "-t " + f"{dic['time_data']}",
        "-w " + f"{dic['dt_data']}",
        "-f " + f"{dic['subfolders']}",
        "-s " + f"{dic['showpywarn']}",
    ]
    print(
        "\nData: Generation of csv files following the SPE11 benchmark format, please wait."
    )
    prosc = subprocess.run(data_exe, check=True)
    if prosc.returncode != 0:
        raise ValueError(f"Invalid result: { prosc.returncode }")
