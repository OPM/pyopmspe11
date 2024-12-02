# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""
Utiliy functions for the simulations, data processing, and plotting.
"""
import os
import subprocess


def simulations(dic, deck, folder, dryrun):
    """
    Run OPM Flow

    Args:
        dic (dict): Global dictionary\n
        deck (str): Name of the input deck\n
        folder (str): Name of destination of the output files\n
        dryrun (bool): Run a dry simulation for the initial files.

    Returns:
        None

    """
    flag = ""
    flow = dic["flow"]
    if dryrun:
        flow = dic["only_flow"]
        flag = "--enable-dry-run=1"
    os.system(
        f"{flow} --output-dir={dic['exe']}/{dic['fol']}/{folder} "
        f"{dic['exe']}/{dic['fol']}/deck/{deck}.DATA {flag} & wait\n"
    )


def plotting(dic):
    """
    Generate the figures

    Args:
        dic (dict): Global dictionary

    Returns:
        None

    """
    os.chdir(f"{dic['exe']}")
    plot_exe = [
        "python3",
        f"{dic['pat']}/visualization/plotting.py",
        "-p " + f"{dic['fol']}",
        "-d " + f"{dic['spe11']}",
        "-g " + f"{dic['generate']}",
        "-r " + f"{dic['resolution']}",
        "-s " + f"{dic['showpywarn']}",
        "-l " + f"{dic['latex']}",
    ]
    print(" ".join(plot_exe))
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
    os.chdir(f"{dic['exe']}")
    data_exe = [
        "python3",
        f"{dic['pat']}/visualization/data.py",
        "-p " + f"{dic['fol']}",
        "-d " + f"{dic['spe11']}",
        "-g " + f"{dic['generate']}",
        "-r " + f"{dic['resolution']}",
        "-t " + f"{dic['time_data']}",
        "-w " + f"{dic['dt_data']}",
        "-u " + f"{dic['use']}",
        "-s " + f"{dic['showpywarn']}",
    ]
    print(" ".join(data_exe))
    prosc = subprocess.run(data_exe, check=True)
    if prosc.returncode != 0:
        raise ValueError(f"Invalid result: { prosc.returncode }")
