# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT
# pylint: disable=R0912, R0915

"""Main script for pyopmspe11"""

import os
import sys
import argparse
import subprocess

from pyopmspe11.utils.inputvalues import process_input, check_deck
from pyopmspe11.utils.runs import simulations, plotting, data
from pyopmspe11.visualization.plotting import plot_results
from pyopmspe11.utils.mapproperties import generate_files


def main(argv: list[str] | None = None) -> None:
    """Main entry point"""
    args = load_parser(argv)

    if args["compare"]:
        print("\nCompare: Generating common plots to compare results, please wait.")
        plot_results({"compare": args["compare"]})
        print(f"\nThe figures have been written to {os.getcwd()}/compare/")
        return

    cfg = process_input(args)
    cfg.deckfol = f"{cfg.fol}/deck" if cfg.subfolders == "1" else cfg.fol
    flowfol = f"{cfg.fol}/flow" if cfg.subfolders == "1" else cfg.fol
    make_dir(cfg.fol)
    os.chdir(cfg.fol)

    if cfg.mode == "all" or "deck" in cfg.mode:
        if cfg.subfolders == "1":
            make_dir(cfg.deckfol)
        print("\nDeck: Generating the input files, please wait.")
        generate_files(cfg)
        print(f"\nThe deck files have been written to {cfg.deckfol}")

    if cfg.mode == "all" or "flow" in cfg.mode:
        check_deck(cfg)
        print("\nFlow: Running the simulations, please wait.")
        simulations(cfg, flowfol)
        print(f"\nThe simulation results have been written to {flowfol}")

    if cfg.mode == "all" or "data" in cfg.mode:
        make_dir(f"{cfg.fol}/data" if cfg.subfolders == "1" else cfg.fol)
        data(cfg)

    if cfg.mode == "all" or "plot" in cfg.mode:
        make_dir(f"{cfg.fol}/figures" if cfg.subfolders == "1" else cfg.fol)
        plotting(cfg)


def make_dir(path: str) -> None:
    """Create directory if missing"""
    if not os.path.exists(path):
        subprocess.run(["mkdir", "-p", path], check=True)


def load_parser(argv: list[str] | None) -> dict:
    """CLI arguments"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="pyopmspe11, a Python tool for the three SPE11 benchmark"
        " cases provided by the Open Porous Media (OPM) project.",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str.strip,
        default="input.toml",
        help="The base name of the input file",
    )
    parser.add_argument(
        "-m",
        "--mode",
        type=str.strip,
        choices=[
            "deck",
            "flow",
            "data",
            "plot",
            "deck_flow",
            "flow_data",
            "data_plot",
            "deck_flow_data",
            "flow_data_plot",
            "all",
        ],
        default="deck_flow",
        help="Parts of pyopmsoe11 to run",
    )
    parser.add_argument(
        "-c",
        "--compare",
        type=str.strip,
        choices=["spe11a", "spe11b", "spe11c", ""],
        default="",
        help="Generate a common plot for the current folders",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str.strip,
        default="output",
        help="The base name of the output folder",
    )
    parser.add_argument(
        "-t",
        "--time",
        type=str.strip,
        default="5",
        help="If one number, time step for the spatial maps (spe11a [h]; spe11b/c "
        "[y]); otherwise, times separated by commas",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=str.strip,
        default="8,1,5",
        help="Number of x, y, and z elements to map the simulation results to the "
        "dense report data",
    )
    parser.add_argument(
        "-g",
        "--generate",
        type=str.strip,
        default="performance_sparse",
        choices=[
            "dense",
            "sparse",
            "performance",
            "performance-spatial",
            "dense_performance",
            "dense_sparse",
            "performance_sparse",
            "dense_performance-spatial",
            "dense_performance_sparse",
            "all",
        ],
        help="Type of fata to generate",
    )
    parser.add_argument(
        "-w",
        "--write",
        type=str.strip,
        default="0.1",
        help="Time interval for the sparse and performance data (spe11a [h]; spe11b/c [y])",
    )
    parser.add_argument(
        "-f",
        "--subfolders",
        choices=["0", "1"],
        type=str.strip,
        default="1",
        help="Set to 0 to not create the subfolders deck, flow, data, and figures, i.e., to "
        "write all generated files in the output directory",
    )
    parser.add_argument(
        "-n",
        "--neighbourhood",
        choices=["lower", ""],
        type=str.strip,
        default="",
        help="Region to model (the default '' means the whole system)",
    )
    return vars(parser.parse_known_args(argv)[0])


if __name__ == "__main__":
    main(sys.argv[1:])
