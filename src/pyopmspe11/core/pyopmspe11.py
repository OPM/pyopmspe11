# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT
# pylint: disable=R0912, R0915

"""Main script for pyopmspe11"""

import argparse
import os
import subprocess
import sys

from pyopmspe11.utils.inputvalues import check_deck, process_input
from pyopmspe11.utils.mapproperties import generate_files
from pyopmspe11.utils.runs import data, plotting, simulations
from pyopmspe11.visualization.plotting import plot_results


def main(argv: list[str] | None = None) -> None:
    """Main entry point"""
    args = load_parser(argv)
    check_cmdargs(args)

    if args.compare:
        print("\nCompare: Generating common plots to compare results, please wait.")
        plot_results({"compare": args.compare})
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


def load_parser(argv: list[str] | None) -> argparse.Namespace:
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
        help="Parts of pyopmspe11 to run",
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
        help="Type of data to generate",
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
    return parser.parse_args(argv)


def check_cmdargs(cmdargs: argparse.Namespace) -> None:
    """Validate command-line arguments and incompatible operations.

    The checks cover configuration and output names, spatial resolution,
    reporting times and intervals, comparison mode, and options restricted
    to data-generation workflows.

    Parameters
    ----------
    cmdargs
        Parsed arguments returned by :func:`load_parser`.

    Raises
    ------
    SystemExit
        If an argument is invalid or an incompatible combination is requested.
    """
    input_file = cmdargs.input
    if not input_file:
        print("\nInvalid value for '-i', the input file cannot be empty.\n")
        raise SystemExit(1)
    if not input_file.lower().endswith((".toml", ".txt")):
        print(
            f"\nInvalid extension for input file '-i {input_file}', "
            "valid extensions are .toml or .txt.\n"
        )
        raise SystemExit(1)
    if not cmdargs.output:
        print("\nInvalid value for '-o', the output folder cannot be empty.\n")
        raise SystemExit(1)
    resolution = cmdargs.resolution
    try:
        resolution_values = [int(value.strip()) for value in resolution.split(",")]
    except ValueError:
        resolution_values = []
    if len(resolution_values) != 3 or any(value <= 0 for value in resolution_values):
        print(
            f"\nInvalid value '-r {resolution}', expected three positive "
            "integers separated by commas, e.g., '-r 8,1,5'.\n"
        )
        raise SystemExit(1)
    time = cmdargs.time
    try:
        time_values = [float(value.strip()) for value in time.split(",")]
    except ValueError:
        time_values = []
    if not time_values or any(value < 0 for value in time_values):
        print(
            f"\nInvalid value '-t {time}', expected non-negative numbers "
            "separated by commas.\n"
        )
        raise SystemExit(1)
    write = cmdargs.write
    try:
        write_value = float(write)
    except ValueError:
        write_value = 0
    if write_value <= 0:
        print(f"\nInvalid value '-w {write}', expected a positive number.\n")
        raise SystemExit(1)
    mode = cmdargs.mode
    has_data = mode == "all" or "data" in mode
    data_options = {
        "-g": ("generate", "performance_sparse"),
        "-r": ("resolution", "8,1,5"),
        "-t": ("time", "5"),
        "-w": ("write", "0.1"),
    }
    if not has_data:
        invalid_options = [
            option
            for option, (name, default) in data_options.items()
            if getattr(cmdargs, name) != default
        ]
        if invalid_options:
            print(
                f"\nInvalid option for '-m {mode}'; {', '.join(invalid_options)} "
                "can only be used when the selected mode writes benchmark "
                "data.\n"
            )
            raise SystemExit(1)
    compare = cmdargs.compare
    if compare:
        compare_options = {
            "-i": ("input", "input.toml"),
            "-m": ("mode", "deck_flow"),
            "-o": ("output", "output"),
            "-t": ("time", "5"),
            "-r": ("resolution", "8,1,5"),
            "-g": ("generate", "performance_sparse"),
            "-w": ("write", "0.1"),
            "-f": ("subfolders", "1"),
            "-n": ("neighbourhood", ""),
        }
        invalid_options = [
            option
            for option, (name, default) in compare_options.items()
            if getattr(cmdargs, name) != default
        ]
        if invalid_options:
            print(
                "\nInvalid combination, '-c' runs the standalone comparison "
                "workflow and cannot be combined with "
                f"{', '.join(invalid_options)}.\n"
            )
            raise SystemExit(1)


if __name__ == "__main__":
    main(sys.argv[1:])
