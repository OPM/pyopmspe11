# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""
Utiliy functions to run the studies.
"""
import os
import subprocess


def simulations(dic, deck, folder):
    """
    Function to run OPM Flow

    Args:
        dic (dict): Global dictionary with required parameters
        deck: Name of the input deck
        folder: destination of the output files

    """
    os.system(
        f"{dic['flow']} --output-dir={dic['exe']}/{dic['fol']}/{folder} "
        f"{dic['exe']}/{dic['fol']}/deck/{deck}.DATA  & wait\n"
    )


def plotting(dic):
    """
    Function to generate and run the plotting.py file

    Args:
        dic (dict): Global dictionary with required parameters

    """
    os.chdir(f"{dic['exe']}")
    plot_exe = [
        "python3",
        f"{dic['pat']}/visualization/plotting.py",
        "-p " + f"{dic['fol']}",
        "-d " + f"{dic['spe11']}",
        "-g " + f"{dic['generate']}",
        "-r " + f"{dic['resolution']}",
        "-t " + f"{dic['time_data']}",
    ]
    print(" ".join(plot_exe))
    prosc = subprocess.run(plot_exe, check=True)
    if prosc.returncode != 0:
        raise ValueError(f"Invalid result: { prosc.returncode }")


def data(dic):
    """
    Function to write the sparse and dense benchmark data

    Args:
        dic (dict): Global dictionary with required parameters

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
        "-l " + f"{dic['load']}",
        "-u " + f"{dic['use']}",
    ]
    print(" ".join(data_exe))
    prosc = subprocess.run(data_exe, check=True)
    if prosc.returncode != 0:
        raise ValueError(f"Invalid result: { prosc.returncode }")
