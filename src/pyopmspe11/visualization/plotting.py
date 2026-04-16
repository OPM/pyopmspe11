# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT
# pylint: disable=R0902,R0912,R0801,R0913,R0914,R0915,R0917

"""Script to plot the results"""

import os
import argparse
import shutil
import math as mt
import subprocess
from io import StringIO
from dataclasses import dataclass
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors, ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

SECONDS_IN_YEAR = 31536000.0

font = {"family": "normal", "weight": "normal", "size": 20}
matplotlib.rc("font", **font)


@dataclass(frozen=True)
class CaseConfig:
    """Case configuration"""

    case: str
    tlabel: str
    dims: int
    tscale: float
    lower: bool


@dataclass(frozen=True)
class RunConfig:
    """Runtime configuration"""

    folders: list
    generate: str
    compare: str
    where: str
    dataf: str
    colors: list
    linestyles: list
    props: dict


@dataclass(frozen=True)
class GridState:
    """Grid configuration"""

    times: list
    xmsh: np.ndarray
    zmsh: np.ndarray
    xmx: list
    ymy: list
    zmz: list
    kinds: list
    cmaps: list
    dims: int


def configure_matplotlib():
    """Parameters for the figures"""
    if shutil.which("latex") == "None":
        print("\nLaTeX is recommended for the figures.")
    plt.rcParams.update(
        {
            "text.usetex": shutil.which("latex") != "None",
            "font.family": "monospace",
            "legend.columnspacing": 0.9,
            "legend.handlelength": 3.5,
            "legend.fontsize": 15,
            "lines.linewidth": 4,
            "axes.titlesize": 20,
            "axes.grid": True,
            "figure.figsize": (10, 5),
        }
    )


def load_csv(path):
    """Load the csv"""
    return np.genfromtxt(path, delimiter=",", skip_header=1)


def load_time_series(folder, case):
    """Read the time_series.csv"""
    path = f"{folder}/data/{case}_time_series.csv"
    if not os.path.isfile(path):
        path = f"{folder}/{case}_time_series.csv"
    return load_csv(path)


def load_performance(folder, case, kind):
    """Read the performance_time_series*.csv"""
    path = f"{folder}/data/{case}_performance_time_series{kind}.csv"
    if not os.path.isfile(path):
        path = f"{folder}/{case}_performance_time_series{kind}.csv"
    return load_csv(path)


def load_spatial(folder, dataf, case, kind, time, tlabel):
    """Read the spatial_map_*.csv"""
    return load_csv(f"{folder}{dataf}/{case}{kind}_spatial_map_{time}{tlabel}.csv")


def performance_label(csv, index, folder):
    """Set the metrics in the labels for the performance plots"""
    stats = [
        f"sum={np.sum(csv[:,1]):.3e}",
        f"sum={np.sum(csv[:,2]):.3e}",
        f"max={np.max(csv[:,3]):.3e}",
        f"max={csv[-1,4]:.3e}",
        f"sum={np.sum(csv[:,5]):.3e}",
        f"sum={np.sum(csv[:,6]):.3e}",
        f"sum={np.sum(csv[:,7]):.3e}",
        f"sum={np.sum(csv[:,8]):.3e}",
        f"sum={np.sum(csv[:,9]):.3e}",
    ]
    return f"{stats[index]} ({folder})"


def generate_grid(folder, dataf, tlabel, dims, time):
    """Create the meshgrid"""
    files = [f for f in os.listdir(f"{folder}{dataf}") if f.endswith(f"{tlabel}.csv")]
    times = np.array([int(f[19:-5]) for f in files if len(f) < 30])
    if times.size == 0:
        times = np.array([int(f[31:-5]) for f in files])
    times = list(times[np.argsort(times)])
    if time.size == 1:
        if time > 0:
            times = list(range(0, times[-1] + 1, time))
        else:
            times = [time]
    else:
        times = list(time)
    csv = load_csv(f"{folder}{dataf}/{files[0]}")
    length = csv[-1][0] + csv[0][0]
    width = csv[-1][dims - 2] + csv[0][dims - 2]
    height = csv[-1][dims - 1] + csv[0][dims - 1]
    xmx = np.linspace(0, length, round(length / (2.0 * csv[0][0])) + 1)
    ymy = np.linspace(0, width, round(width / (2.0 * csv[0][dims - 2])) + 1)
    zmz = np.linspace(0, height, round(height / (2.0 * csv[0][dims - 1])) + 1)
    xmsh, zmsh = np.meshgrid(xmx, zmz[::-1])
    return times, xmsh, zmsh, xmx, ymy, zmz


def performance(case_cfg, run_cfg):
    """Plot the performance"""
    for kind in ["", "_detailed"]:
        fig = plt.figure(figsize=(40, 75))
        plots = [
            "tstep",
            "fsteps",
            "mass",
            "dof",
            "nliter",
            "nres",
            "liniter",
            "runtime",
            "tlinsol",
        ]
        ylabels = ["s", "\\#", "kg", "\\#", "\\#", "\\#", "\\#", "s", "s"]
        for i, (plot, ylabel) in enumerate(zip(plots, ylabels)):
            axis = fig.add_subplot(9, 5, i + 1)
            for folder_index, folder in enumerate(run_cfg.folders):
                csv = load_performance(folder, case_cfg.case, kind)
                if len(csv.flatten()) < 12:
                    csv = np.array([csv])
                times = csv[:, 0] / case_cfg.tscale
                label = performance_label(csv, i, folder.split("/")[-1])
                axis.step(
                    times,
                    csv[:, i + 1],
                    lw=2,
                    color=run_cfg.colors[folder_index],
                    label=label,
                )
            axis.set_title(f"{plot}, {case_cfg.case}")
            axis.set_ylabel(ylabel)
            axis.set_xlabel(f"Time [{case_cfg.tlabel}]")
            axis.legend()
        fig.savefig(
            f"{run_cfg.where}/{case_cfg.case}_performance{kind}.png",
            bbox_inches="tight",
        )


def sparse_data(case_cfg, run_cfg):
    """Plot the sparse data"""
    fig = plt.figure(figsize=(25, 40))
    plots = ["sensors", "boxA", "boxB", "boxC", "facie 1"]
    ylabels = ["Presure [Pa]", "Mass [kg]", "Mass [kg]", "Length [m]", "Mass [kg]"]
    labels = [
        ["p1", "p2"],
        ["mobA", "immA", "dissA", "sealA"],
        ["mobB", "immB", "dissB", "sealB"],
        ["MC"],
        ["sealTot"],
    ]
    nfigs = 5
    if case_cfg.case != "spe11a":
        plots.append("boundaries")
        ylabels.append("Mass [kg]")
        labels.append(["boundTot"])
        nfigs += 1
    for k, (plot, ylabel) in enumerate(zip(plots, ylabels)):
        axis = fig.add_subplot(nfigs, 3, k + 1)
        for folder_index, folder in enumerate(run_cfg.folders):
            column = sum(len(labels[i]) for i in range(k)) + 1
            axis.text(
                0.7,
                0.15 + folder_index * 0.05,
                run_cfg.folders[-1 - folder_index].split("/")[-1],
                transform=axis.transAxes,
                verticalalignment="top",
                bbox=run_cfg.props,
                color=run_cfg.colors[len(run_cfg.folders) - folder_index - 1],
            )
            csv = load_time_series(folder, case_cfg.case)
            times = csv[:, 0] / case_cfg.tscale
            for j, label in enumerate(labels[k]):
                if folder_index == 0:
                    axis.plot(
                        times,
                        csv[:, column],
                        label=label,
                        color="k",
                        linestyle=run_cfg.linestyles[-1 + j],
                    )
                if label == "MC":
                    axis.plot(
                        times,
                        csv[:, column],
                        color=run_cfg.colors[folder_index],
                        linestyle=run_cfg.linestyles[-1 + j + folder_index],
                    )
                else:
                    axis.plot(
                        times,
                        csv[:, column],
                        color=run_cfg.colors[folder_index],
                        linestyle=run_cfg.linestyles[-1 + j],
                    )
                column += 1
        axis.set_title(f"{plot}, {case_cfg.case}")
        axis.set_ylabel(ylabel)
        axis.set_xlabel(f"Time [{case_cfg.tlabel}]")
        axis.legend()
    fig.savefig(f"{run_cfg.where}/{case_cfg.case}_sparse_data.png", bbox_inches="tight")


def dense_data(case_cfg, run_cfg, grid):
    """Plot the dense data"""
    for kind in grid.kinds:
        if kind == "":
            quantities = [
                "pressure",
                "sgas",
                "xco2",
                "xh20",
                "gden",
                "wden",
                "tco2",
                "temp",
            ]
            units = [
                "[Pa]",
                "[-]",
                "[-]",
                "[-]",
                r"[kg/m$^3$]",
                r"[kg/m$^3$]",
                "[kg]",
                "C",
            ]
            allplots = [-1] * 8
        else:
            quantities = [
                "cvol",
                "arat",
                "CO2 max_norm_res",
                "H2O max_norm_res",
                "CO2 mb_error",
                "H2O mb_error",
            ]
            units = [r"[m$^3$]", "[-]", "[-]", "[-]", "[-]", "[-]"]
            allplots = [0, 0, -1, -1, -1, -1]
        for qi, quantity in enumerate(quantities):
            csvs = [
                load_spatial(
                    run_cfg.folders[0],
                    run_cfg.dataf,
                    case_cfg.case,
                    kind,
                    t,
                    case_cfg.tlabel,
                )
                for t in grid.times
            ]
            if qi == csvs[0].shape[1] - grid.dims:
                break
            first = csvs[0][:, grid.dims + qi]
            if np.isnan(first).all():
                continue
            minc, maxc = np.nanmin(first), np.nanmax(first)
            ptimes = grid.times[: allplots[qi]] + [grid.times[-1]]
            if case_cfg.case != "spe11a":
                fig = plt.figure(figsize=(50, 3 * len(ptimes)))
                if case_cfg.lower:
                    fig = plt.figure(figsize=(100, 3 * len(ptimes)))
            else:
                fig = plt.figure(figsize=(45, 6.5 * len(ptimes)), dpi=80)
            plots = []
            mins = []
            maxs = []
            sums = []
            for ti, time in enumerate(ptimes):
                values = csvs[ti][:, grid.dims + qi]
                mins.append(np.nanmin(values))
                maxs.append(np.nanmax(values))
                if quantity == "tco2":
                    sums.append(np.sum(values[values >= 0]))
                minc = min(minc, mins[-1])
                maxc = max(maxc, maxs[-1])
                arr = np.zeros((len(grid.zmz) - 1, len(grid.xmx) - 1))
                for zi in range(len(grid.zmz) - 1):
                    if case_cfg.case != "spe11c":
                        start = zi * (len(grid.xmx) - 1)
                        arr[-1 - zi, :] = values[start : start + (len(grid.xmx) - 1)]
                    else:
                        mid = mt.floor((len(grid.ymy) - 1) / 2)
                        slice_size = len(grid.xmx) - 1
                        offset = (zi * (len(grid.ymy) - 1) + mid) * slice_size
                        arr[-1 - zi, :] = values[offset : offset + slice_size]
                plots.append(arr)
            for j, time in enumerate(ptimes):
                axis = fig.add_subplot(len(ptimes), 3, j + 1)
                imag = axis.pcolormesh(
                    grid.xmsh,
                    grid.zmsh,
                    plots[j],
                    shading="flat",
                    cmap=(
                        grid.cmaps[qi]
                        if mins != maxs
                        else colors.ListedColormap(["#1319bf"])
                    ),
                )
                if quantity == "tco2":
                    title = f"{time}{case_cfg.tlabel}, {quantity} {units[qi]}(sum={sums[j]:.1E})"
                else:
                    prefix = f"{time}{case_cfg.tlabel}, " if allplots[qi] == -1 else ""
                    title = f"{prefix}{quantity} {units[qi]}(min={mins[j]:.1E}, max={maxs[j]:.1E})"
                axis.set_title(
                    f"{title}, {case_cfg.case} ({run_cfg.folders[0].split('/')[-1]})"
                )
                axis.axis("scaled")
                axis.xaxis.set_major_locator(ticker.MaxNLocator(14))
                axis.yaxis.set_major_locator(ticker.MaxNLocator(4))
                imag.set_clim(minc, maxc)
                if j % 3 != 0:
                    axis.set_yticks([])
                if (
                    j
                    < ((len(ptimes) - len(ptimes) % 3) / 3) * 3 - (3 - len(ptimes) % 3)
                    or (len(ptimes) % 3 == 1 and j == len(ptimes) - 4)
                    or (len(ptimes) % 3 == 2 and j == len(ptimes) - 5)
                ):
                    axis.set_xticks([])
                if (
                    (j + 1) % 3 == 0
                    or len(ptimes) == 1
                    or (len(ptimes) == 2 and j == 1)
                ):
                    divider = make_axes_locatable(axis)
                    fig.colorbar(
                        imag,
                        cax=divider.append_axes("right", size="5%", pad=0.05),
                        ticks=np.linspace(minc, maxc, 5),
                        format=lambda x, _: f"{x:.2e}",
                    )
                if case_cfg.lower:
                    axis.set_ylim([0, 0.55] if case_cfg.case == "spe11a" else [0, 550])
            fig.savefig(
                f"{run_cfg.where}/{case_cfg.case}_{quantity}_2Dmaps.png",
                bbox_inches="tight",
            )
            plt.close(fig)


def plot_results(args):
    """Orchestrate the plotting"""
    configure_matplotlib()
    where = ""
    dataf = ""
    if args["compare"]:
        args["deck"] = args["compare"]
        args["neighbourhood"] = ""
        where = "compare/"
        folders = sorted(
            [n for n in os.listdir(".") if os.path.isdir(n) and n != "compare"]
        )
        if not os.path.isdir("compare"):
            subprocess.run(["mkdir", "compare"], check=False)
    else:
        folders = [args["folder"].strip()]
        if int(args["subfolders"]) == 1:
            dataf = "/data"
            where = f"{folders[0]}/figures"
        else:
            where = folders[0]
    if args["deck"] == "spe11a":
        case_cfg = CaseConfig(args["deck"], "h", 2, 3600.0, bool(args["neighbourhood"]))
    else:
        case_cfg = CaseConfig(
            args["deck"], "y", 2, SECONDS_IN_YEAR, bool(args["neighbourhood"])
        )
    if args["deck"] == "spe11c":
        case_cfg = CaseConfig(
            args["deck"], "y", 3, SECONDS_IN_YEAR, bool(args["neighbourhood"])
        )
    run_cfg = RunConfig(
        folders,
        args["generate"],
        args["compare"],
        where,
        dataf,
        [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
            "r",
            "k",
        ],
        [
            "--",
            (0, (1, 1)),
            "-.",
            (0, (1, 10)),
            (0, (1, 1)),
            (5, (10, 3)),
            (0, (5, 10)),
            (0, (5, 5)),
            (0, (5, 1)),
            (0, (3, 10, 1, 10)),
            (0, (3, 5, 1, 5)),
            (0, (3, 1, 1, 1)),
            (0, (3, 5, 1, 5, 1, 5)),
            (0, (3, 10, 1, 10, 1, 10)),
            (0, (3, 1, 1, 1, 1, 1)),
            (0, ()),
            "-",
        ],
        {"boxstyle": "round", "facecolor": "wheat", "alpha": 0.1},
    )
    if args["generate"] in [
        "all",
        "performance",
        "dense_performance",
        "performance_sparse",
        "dense_performance_sparse",
    ]:
        performance(case_cfg, run_cfg)
    if args["generate"] in [
        "all",
        "sparse",
        "dense_sparse",
        "performance_sparse",
        "dense_performance_sparse",
    ]:
        sparse_data(case_cfg, run_cfg)
    if args["compare"]:
        return
    plt.rcParams.update({"axes.grid": False})
    if args["generate"] in [
        "all",
        "dense",
        "performance-spatial",
        "dense_performance",
        "dense_sparse",
        "dense_performance-spatial",
        "dense_performance_sparse",
    ]:
        time = np.genfromtxt(StringIO(args["time"]), delimiter=",", dtype=int)
        times, xmsh, zmsh, xmx, ymy, zmz = generate_grid(
            run_cfg.folders[0], run_cfg.dataf, case_cfg.tlabel, case_cfg.dims, time
        )
        kinds = (
            ["", "_performance"]
            if args["generate"] in ["all", "dense_performance-spatial"]
            else [""] if args["generate"].startswith("dense") else ["_performance"]
        )
        grid = GridState(
            times,
            xmsh,
            zmsh,
            xmx,
            ymy,
            zmz,
            kinds,
            [
                "seismic",
                "jet",
                "viridis",
                "viridis_r",
                "PuOr",
                "PuOr_r",
                "turbo",
                "coolwarm",
            ],
            case_cfg.dims,
        )
        dense_data(case_cfg, run_cfg, grid)


def main():
    """Entry point"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Main script to plot the results",
    )
    parser.add_argument("-p", "--folder", default="output", type=str.strip)
    parser.add_argument("-c", "--compare", default="", type=str.strip)
    parser.add_argument("-d", "--deck", default="spe11b", type=str.strip)
    parser.add_argument("-g", "--generate", default="sparse", type=str.strip)
    parser.add_argument("-f", "--subfolders", default="1", type=str.strip)
    parser.add_argument("-t", "--time", default="5", type=str.strip)
    parser.add_argument("-n", "--neighbourhood", default="", type=str.strip)
    args = vars(parser.parse_known_args()[0])
    plot_results(args)
    print(f"The png figures have been saved on {args['folder']}")


if __name__ == "__main__":
    main()
