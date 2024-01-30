# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Main script"""
import os
import argparse
from pyopmspe11.utils.inputvalues import process_input, check_deck
from pyopmspe11.utils.runs import simulations, plotting, data
from pyopmspe11.visualization.plotting import plot_results
from pyopmspe11.utils.writefile import opm_files, initial
from pyopmspe11.utils.mapproperties import (
    grid,
    positions,
)


def pyopmspe11():
    """Main function"""
    cmdargs = load_parser()
    file = cmdargs["input"].strip()  # Name of the input file
    dic = {"fol": cmdargs["output"].strip()}  # Name for the output folder
    dic["generate"] = cmdargs["generate"].strip()  # What data to write
    dic["exe"] = os.getcwd()  # Path to the folder of the input.txt file
    dic["mode"] = cmdargs["mode"].strip()  # Parts of the workflow to run
    dic["pat"] = os.path.dirname(__file__)[:-5]  # Path to the pyopmspe11 folder
    dic["compare"] = cmdargs["compare"].strip()  # Make common figures for comparison
    dic["use"] = cmdargs["use"].strip()  # OPM or resdata python package
    dic["resolution"] = cmdargs[
        "resolution"
    ].strip()  # Spatial resolution to write the data
    dic["time_data"] = cmdargs["time"]  # Temporal resolution to write the dense data
    dic["dt_data"] = float(
        cmdargs["write"].strip()
    )  # Temporal resolution to write the sparse and performance data
    # If the compare plots are generated, then we exit right afterwards
    if dic["compare"]:
        plot_results(dic)
        return

    # Process the input file (open pyopmspe11.utils.inputvalues to see the abbreviations meaning)
    dic = process_input(dic, file)

    # Make the output folders
    if not os.path.exists(f"{dic['exe']}/{dic['fol']}"):
        os.system(f"mkdir {dic['exe']}/{dic['fol']}")
    for fil in ["deck", "flow"]:
        if not os.path.exists(f"{dic['exe']}/{dic['fol']}/{fil}"):
            os.system(f"mkdir {dic['exe']}/{dic['fol']}/{fil}")
    os.chdir(f"{dic['exe']}/{dic['fol']}")

    if dic["mode"] in ["all", "deck", "deck_flow", "deck_flow_data"]:
        # Initialize the grid
        dic = grid(dic)

        # For corner-point grids, get the cell centers by executing flow
        if dic["grid"] == "corner-point":
            initial(dic)
            os.chdir(f"{dic['exe']}/{dic['fol']}/deck")
            simulations(dic, "INITIAL", "flow")
            print("Files used to generate the corner-point grid (INITIAL.* files)")

        # Check the generated deck, flow version, and chosen co2store implementation
        check_deck(dic)

        # Get the sand and well positions
        dic = positions(dic)

        # Write used opm related files
        opm_files(dic)

    if dic["mode"] in ["all", "flow", "deck_flow", "deck_flow_data"]:
        # Run the simulations
        simulations(dic, dic["fol"].upper(), "flow")

    if dic["mode"] in ["all", "data", "deck_flow_data", "data_plot"]:
        # Write the data
        if not os.path.exists(f"{dic['exe']}/{dic['fol']}/data"):
            os.system(f"mkdir {dic['exe']}/{dic['fol']}/data")
        data(dic)

    if dic["mode"] in ["all", "plot", "data_plot"]:
        # Make some useful plots after the studies
        if not os.path.exists(f"{dic['exe']}/{dic['fol']}/figures"):
            os.system(f"mkdir {dic['exe']}/{dic['fol']}/figures")
        plotting(dic)


def load_parser():
    """Argument options"""
    parser = argparse.ArgumentParser(
        description="Main script to run the spe11s with OPM Flow."
    )
    parser.add_argument(
        "-i",
        "--input",
        default="input.txt",
        help="The base name of the input file ('input.txt' by default).",
    )
    parser.add_argument(
        "-m",
        "--mode",
        default="deck_flow",
        help="Run the whole framework ('all'), only create decks ('deck'), "
        "only run flow ('flow'), only write benchmark data ('data'), "
        "only create plots ('plot'), deck and run ('deck_flow'), "
        "data and plot (data_plot), or deck, run, and data "
        "(deck_flow_data) ('deck_flow' by default).",
    )
    parser.add_argument(
        "-c",
        "--compare",
        default="",
        help="Generate a common plot for the current folders ('' by default).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="output",
        help="The base name of the output folder ('output' by default).",
    )
    parser.add_argument(
        "-t",
        "--time",
        default="5",
        help="If one number, time step for the spatial maps (spe11a [h]; spe11b/c "
        "[y]) ('5' by default); otherwise, times separated by commas.",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        default="8,1,5",
        help="Number of x, y, and z elements to write the data ('8,1,5' by default).",
    )
    parser.add_argument(
        "-g",
        "--generate",
        default="performance_sparse",
        help="Write only the 'dense', 'sparse', 'performance', 'dense_performance', "
        "'performance_sparse', 'dense_sparse', or 'all'",
    )
    parser.add_argument(
        "-u",
        "--use",
        default="resdata",
        help="Using the 'opm' or 'resdata' python package (resdata by default).",
    )
    parser.add_argument(
        "-w",
        "--write",
        default="0.1",
        help="Time interval for the sparse and performance data (spe11a [h]; spe11b/c [y]) "
        "('0.1' by default).",
    )
    return vars(parser.parse_known_args()[0])


def main():
    """Main function"""
    pyopmspe11()
