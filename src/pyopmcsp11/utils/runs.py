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
        f"{dic['exe']}/{dic['fol']}/preprocessing/{deck}.DATA  & wait\n"
    )


def plotting(dic, time):
    """
    Function to generate and run the plotting.py file

    Args:
        dic (dict): Global dictionary with required parameters

    """
    os.chdir(f"{dic['exe']}/{dic['fol']}/postprocessing")
    prosc = subprocess.run(
        [
            "python",
            f"{dic['pat']}/visualization/plotting.py",
            f"-t {time}",
            "-p " + f"{dic['exe']}/{dic['fol']}",
            "-d " + f"{dic['csp11'].upper()}",
        ],
        check=True,
    )
    if prosc.returncode != 0:
        raise ValueError(f"Invalid result: { prosc.returncode }")
