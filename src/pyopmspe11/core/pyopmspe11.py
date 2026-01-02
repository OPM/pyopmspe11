# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT
# pylint: disable=R0912, R0915

"""Main script for pyopmspe11"""

import os
import argparse
import warnings
from pyopmspe11.utils.inputvalues import process_input, check_deck
from pyopmspe11.utils.runs import simulations, plotting, data
from pyopmspe11.visualization.plotting import plot_results
from pyopmspe11.utils.writefile import opm_files
from pyopmspe11.utils.mapproperties import grid, positions


def pyopmspe11():
    """Main function for the pyopmspe11 executable"""
    cmdargs = load_parser()
    file = cmdargs["input"].strip()  # Name of the input file
    dic = {"fol": os.path.abspath(cmdargs["output"])}  # Name for the output folder
    dic["generate"] = cmdargs["generate"].strip()  # What data to write
    dic["mode"] = cmdargs["mode"].strip()  # Parts of the workflow to run
    dic["pat"] = os.path.dirname(__file__)[:-5]  # Path to the pyopmspe11 folder
    dic["compare"] = cmdargs["compare"].strip()  # Make common figures for comparison
    dic["resolution"] = cmdargs[
        "resolution"
    ].strip()  # Spatial resolution to write the data
    dic["time_data"] = cmdargs["time"]  # Temporal resolution to write the dense data
    dic["dt_data"] = float(
        cmdargs["write"].strip()
    )  # Temporal resolution to write the sparse and performance data
    dic["showpywarn"] = int(cmdargs["showpywarn"])  # Show or hidde python warnings
    dic["lower"] = bool(cmdargs["neighbourhood"].strip())  # Lower model
    dic["subfolders"] = int(cmdargs["subfolders"])  # Create subfolders
    if dic["showpywarn"] != 1:
        warnings.warn = lambda *args, **kwargs: None
    # If the compare plots are generated, then we exit right afterwards
    if dic["compare"]:
        print("\nCompare: Generating common plots to compare results, please wait.")
        plot_results(dic)
        print(f"\nThe figures have been written to {os.getcwd()}/compare/")
        return

    # Process the input file (open pyopmspe11.utils.inputvalues to see the abbreviations meaning)
    process_input(dic, file)

    # Make the output folders
    if not os.path.exists(f"{dic['fol']}"):
        os.system(f"mkdir {dic['fol']}")
    if dic["subfolders"] == 1:
        dic["deckf"] = f"{dic['fol']}/deck"
        dic["flowf"] = f"{dic['fol']}/flow"
        dic["dataf"] = f"{dic['fol']}/data"
        dic["figsf"] = f"{dic['fol']}/figures"
        if not os.path.exists(dic["deckf"]):
            os.system(f"mkdir {dic['deckf']}")
    else:
        dic["deckf"] = f"{dic['fol']}"
        dic["flowf"] = f"{dic['fol']}"
        dic["dataf"] = f"{dic['fol']}"
        dic["figsf"] = f"{dic['fol']}"
    os.chdir(f"{dic['fol']}")

    if dic["mode"] == "all" or "deck" in dic["mode"]:
        print("\nDeck: Generating the input files, please wait.")
        # Initialize the grid
        grid(dic)
        # Get the sand and well/sources positions
        positions(dic)
        # Write used opm related files
        opm_files(dic)
        print(f"\nThe deck files have been written to {dic['deckf']}")
    if dic["mode"] == "all" or "flow" in dic["mode"]:
        # Check the generated deck and flow version
        check_deck(dic)
        # Run the simulations
        print("\nFlow: Running the simulations, please wait.")
        simulations(dic)
        print(f"\nThe simulation results have been written to {dic['flowf']}")

    if dic["mode"] == "all" or "data" in dic["mode"]:
        # Write the data
        if not os.path.exists(dic["dataf"]):
            os.system(f"mkdir {dic['dataf']}")
        data(dic)

    if dic["mode"] == "all" or "plot" in dic["mode"]:
        # Make some useful plots after the studies
        if not os.path.exists(dic["figsf"]):
            os.system(f"mkdir {dic['figsf']}")
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
        default="input.toml",
        help="The base name of the input file ('input.toml' by default).",
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
        "-f",
        "--subfolders",
        default=1,
        help="Set to 0 to not create the subfolders deck, flow, data, and figures, i.e., to "
        "write all generated files in the output directory ('1' by default).",
    )
    parser.add_argument(
        "-n",
        "--neighbourhood",
        default="",
        help="Region to model; valid options are 'lower' or '' (all reservoir) ('' by default)",
    )
    return vars(parser.parse_known_args()[0])


def main():
    """Main function"""
    pyopmspe11()
