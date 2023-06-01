# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Main script"""
import os
import time
import argparse
from pyopmcsp11.utils.inputvalues import process_input
from pyopmcsp11.utils.runs import simulations, plotting
from pyopmcsp11.utils.writefile import opm_files, initial
from pyopmcsp11.utils.mapproperties import (
    grid,
    positions,
)


def pyopmcsp11():
    """Main function"""
    start_time = time.monotonic()
    parser = argparse.ArgumentParser(
        description="Main script to run the csp11s with OPM Flow."
    )
    parser.add_argument(
        "-i",
        "--input",
        default="input.txt",
        help="The base name of the input file ('input.txt' by default).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="output",
        help="The base name of the output folder ('output' by default).",
    )
    cmdargs = vars(parser.parse_known_args()[0])
    file = cmdargs["input"]  # Name of the input file
    dic = {"fol": cmdargs["output"]}  # Name for the output folder
    dic["pat"] = os.path.dirname(__file__)[:-5]  # Path to the pyopmcsp11 folder
    dic["exe"] = os.getcwd()  # Path to the folder of the input.txt file

    # Process the input file (open pyopmcsp11.utils.inputvalues to see the abbreviations meaning)
    dic = process_input(dic, file)

    # Make the output folders
    if not os.path.exists(f"{dic['exe']}/{dic['fol']}"):
        os.system(f"mkdir {dic['exe']}/{dic['fol']}")
    for fil in ["preprocessing", "jobs", "output"]:
        if not os.path.exists(f"{dic['exe']}/{dic['fol']}/{fil}"):
            os.system(f"mkdir {dic['exe']}/{dic['fol']}/{fil}")
    os.chdir(f"{dic['exe']}/{dic['fol']}")

    # Initialize the grid
    dic = grid(dic)

    # For corner-point grids, get the cell centers by executing flow
    if dic["grid"] == "corner-point":
        initial(dic)
        os.chdir(f"{dic['exe']}/{dic['fol']}/preprocessing")
        simulations(dic, "INITIAL", "preprocessing")

    # Get the sand and well positions
    dic = positions(dic)

    # Write used opm related files
    opm_files(dic)

    # Run the simulations
    simulations(dic, dic["csp11"].upper(), "output")

    # Make some useful plots after the studies (currently only for cartesian grids)
    if dic["grid"] == "cartesian":
        if not os.path.exists(f"{dic['exe']}/{dic['fol']}/postprocessing"):
            os.system(f"mkdir {dic['exe']}/{dic['fol']}/postprocessing")
        plotting(dic, time.monotonic() - start_time)


def main():
    """Main function"""
    pyopmcsp11()
