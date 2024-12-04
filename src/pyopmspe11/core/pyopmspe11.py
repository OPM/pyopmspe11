# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Main script for pyopmspe11"""
import os
import argparse
import warnings
from pyopmspe11.utils.inputvalues import process_input, check_deck, handle_tuning
from pyopmspe11.utils.runs import simulations, plotting, data
from pyopmspe11.visualization.plotting import plot_results
from pyopmspe11.utils.writefile import opm_files, initial
from pyopmspe11.utils.mapproperties import grid, positions


def pyopmspe11():
    """Main function for the pyopmspe11 executable"""
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
    dic["showpywarn"] = int(cmdargs["showpywarn"])  # Show or hidde python warnings
    dic["latex"] = int(cmdargs["latex"])  # LaTeX formatting
    if dic["showpywarn"] != 1:
        warnings.warn = lambda *args, **kwargs: None
    # If the compare plots are generated, then we exit right afterwards
    if dic["compare"]:
        plot_results(dic)
        return

    # Process the input file (open pyopmspe11.utils.inputvalues to see the abbreviations meaning)
    process_input(dic, file)

    # Make the output folders
    if not os.path.exists(f"{dic['exe']}/{dic['fol']}"):
        os.system(f"mkdir {dic['exe']}/{dic['fol']}")
    for fil in ["deck", "flow" if dic["mode"] != "deck" else ""]:
        if not os.path.exists(f"{dic['exe']}/{dic['fol']}/{fil}"):
            os.system(f"mkdir {dic['exe']}/{dic['fol']}/{fil}")
    os.chdir(f"{dic['exe']}/{dic['fol']}")

    if dic["mode"] == "all" or "deck" in dic["mode"]:
        # Check the generated deck, flow version, and chosen co2store implementation
        check_deck(dic)
        # Initialize the grid
        grid(dic)
        # For corner-point grids, get the cell centers by executing flow
        if dic["grid"] == "corner-point":
            initial(dic)
            os.chdir(f"{dic['exe']}/{dic['fol']}/deck")
            simulations(dic, "INITIAL", "flow", True)
            print(
                "Files used to generate the corner-point grid (INITIAL.* files).\n"
                + "Please wait while pyopmspe11 is processing the deck files."
            )
        # Handle tuning
        handle_tuning(dic)
        # Get the sand and well/sources positions
        positions(dic)
        # Write used opm related files
        opm_files(dic)
        print(f"The deck files have been written to {dic['exe']}/{dic['fol']}/deck.")
    if dic["mode"] == "all" or "flow" in dic["mode"]:
        # Run the simulations
        simulations(dic, dic["fol"].upper(), "flow", False)

    if dic["mode"] == "all" or "data" in dic["mode"]:
        # Write the data
        if not os.path.exists(f"{dic['exe']}/{dic['fol']}/data"):
            os.system(f"mkdir {dic['exe']}/{dic['fol']}/data")
        data(dic)

    if dic["mode"] == "all" or "plot" in dic["mode"]:
        # Make some useful plots after the studies
        if not os.path.exists(f"{dic['exe']}/{dic['fol']}/figures"):
            os.system(f"mkdir {dic['exe']}/{dic['fol']}/figures")
        plotting(dic)


def load_parser():
    """Argument options"""
    parser = argparse.ArgumentParser(
        description="pyopmspe11, a Python tool for the three SPE11 benchmark"
        " cases provided by the Open Porous Media (OPM) project.",
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
        "data and plot ('data_plot'), run and data ('flow_data'), deck, "
        "run, and data ('deck_flow_data'), or flow, data, and plot "
        "('flow_data_plot') ('deck_flow' by default).",
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
        help="Number of x, y, and z elements to map the simulation results to the "
        "dense report data ('8,1,5' by default).",
    )
    parser.add_argument(
        "-g",
        "--generate",
        default="performance_sparse",
        help="Write only the 'dense', 'sparse', 'performance', 'performance-spatial', "
        "'dense_performance', 'dense_sparse', 'performance_sparse', "
        "'dense_performance-spatial', 'dense_performance_sparse', or 'all' "
        "('performance_sparse') by default",
    )
    parser.add_argument(
        "-u",
        "--use",
        default="resdata",
        help="Using the 'opm' or 'resdata' python package ('resdata' by default).",
    )
    parser.add_argument(
        "-w",
        "--write",
        default="0.1",
        help="Time interval for the sparse and performance data (spe11a [h]; spe11b/c [y]) "
        "('0.1' by default).",
    )
    parser.add_argument(
        "-s",
        "--showpywarn",
        default=0,
        help="Set to 1 to show Python warnings ('0' by default).",
    )
    parser.add_argument(
        "-l",
        "--latex",
        default=1,
        help="Set to 0 to not use LaTeX formatting ('1' by default).",
    )
    return vars(parser.parse_known_args()[0])


def main():
    """Main function"""
    pyopmspe11()
