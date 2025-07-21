# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT
# pylint: disable=R0912, R0801, R0914, R0915

"""
Script to plot the results
"""

import os
import argparse
import shutil
import warnings
import math as mt
from io import StringIO
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

font = {"family": "normal", "weight": "normal", "size": 20}
matplotlib.rc("font", **font)

SECONDS_IN_YEAR = 31536000.0


def main():
    """Generate figures"""
    parser = argparse.ArgumentParser(description="Main script to plot the results")
    parser.add_argument(
        "-p",
        "--folder",
        default="output",
        help="The folder to the opm simulations.",
    )
    parser.add_argument(
        "-c",
        "--compare",
        default="",
        help="Generate a common plot for the current folders.",
    )
    parser.add_argument(
        "-d",
        "--deck",
        default="spe11b",
        help="The simulated spe case.",
    )
    parser.add_argument(
        "-g",
        "--generate",
        default="sparse",
        help="Plot only the 'dense', 'sparse', 'performance', 'performance-spatial', "
        "'dense_performance', 'dense_sparse', 'performance_sparse', "
        "'dense_performance-spatial', 'dense_performance_sparse', or 'all'",
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
        "-t",
        "--time",
        default="5",
        help="If one number, time step for the spatial maps (spe11a [h]; spe11b/c "
        "[y]) ('5' by default); otherwise, times separated by commas.",
    )
    parser.add_argument(
        "-n",
        "--neighbourhood",
        default="",
        help="Region to model; valid options are 'lower' or '' (all reservoir) ('' by default)",
    )
    cmdargs = vars(parser.parse_known_args()[0])
    if int(cmdargs["showpywarn"]) != 1:  # Show or hidde python warnings
        warnings.warn = lambda *args, **kwargs: None
    dic = {"folders": [cmdargs["folder"].strip()]}
    dic["case"] = cmdargs["deck"].strip()
    dic["generate"] = cmdargs["generate"].strip()
    dic["compare"] = cmdargs["compare"]  # No empty, then the create compare folder
    dic["subfolders"] = int(cmdargs["subfolders"]) == 1  # Create subfolders
    dic["time"] = np.genfromtxt(StringIO(cmdargs["time"]), delimiter=",", dtype=int)
    dic["lower"] = bool(cmdargs["neighbourhood"].strip())  # Lower model
    plot_results(dic)
    print(f"The png figures have been saved on {dic['where']}")


def plot_results(dic):
    """
    Make some figures using the csv files from the benchmark data

    Args:
        dic (dict): Global dictionary

    Returns:
        None

    """
    if shutil.which("latex") == "None":
        print(
            "\nLaTeX is recommended for the figures to show the "
            "nice fonts and given formats. You can install it by "
            "following the instructions in the pyopmspe11's "
            "documentation."
        )
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
    dic["colors"] = [
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
    ]
    dic["linestyle"] = [
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
    ]
    dic["props"] = {"boxstyle": "round", "facecolor": "wheat", "alpha": 0.1}
    if dic["compare"]:
        dic["case"] = dic["compare"]
        dic["where"] = "compare/"
        dic["folders"] = sorted(
            [name for name in os.listdir(".") if os.path.isdir(name)]
        )
        if "compare" not in dic["folders"]:
            os.system("mkdir compare")
        else:
            dic["folders"].remove("compare")
    else:
        if dic["subfolders"]:
            dic["dataf"] = "/data"
            dic["where"] = f"{dic['folders'][0]}/figures"
        else:
            dic["dataf"] = ""
            dic["where"] = dic["folders"][0]
    if dic["case"] == "spe11a":
        dic["tlabel"] = "h"
        dic["dims"] = 2
        dic["tscale"] = 3600.0
    else:
        dic["tlabel"] = "y"
        dic["dims"] = 2
        dic["tscale"] = SECONDS_IN_YEAR
    if dic["case"] == "spe11c":
        dic["dims"] = 3
    if dic["generate"] in [
        "all",
        "performance",
        "dense_performance",
        "performance_sparse",
        "dense_performance_sparse",
    ]:
        performance(dic)
    if dic["generate"] in [
        "all",
        "sparse",
        "dense_sparse",
        "performance_sparse",
        "dense_performance_sparse",
    ]:
        sparse_data(dic)
    if dic["compare"]:
        return
    plt.rcParams.update({"axes.grid": False})
    if dic["generate"] in [
        "all",
        "dense",
        "performance-spatial",
        "dense_performance",
        "dense_sparse",
        "dense_performance-spatial",
        "dense_performance_sparse",
    ]:
        dense_data(dic)


def performance(dic):
    """
    Make the plots related to the performance data (e.g., number of Newton iterations)

    Args:
        dic (dict): Global dictionary

    Returns:
        None

    """
    for kind in ["", "_detailed"]:
        dic["fig"] = plt.figure(figsize=(40, 75))
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
        for k, (plot, ylabel) in enumerate(zip(plots, ylabels)):
            axis = dic["fig"].add_subplot(9, 5, k + 1)
            for nfol, fol in enumerate(dic["folders"]):
                if os.path.isfile(
                    f"{fol}/data/{dic['case']}_performance_time_series{kind}.csv"
                ):
                    csv = np.genfromtxt(
                        f"{fol}/data/{dic['case']}_performance_time_series{kind}.csv",
                        delimiter=",",
                        skip_header=1,
                    )
                else:
                    csv = np.genfromtxt(
                        f"{fol}/{dic['case']}_performance_time_series{kind}.csv",
                        delimiter=",",
                        skip_header=1,
                    )
                labels = [
                    f"sum={np.sum(csv[:,1]):.3e}",
                    f"sum={np.sum(csv[:,2]):.3e}",
                    f"max={np.max(csv[:,3]):.3e}",
                    f"max={csv[-1][4]:.3e}",
                    f"sum={np.sum(csv[:,5]):.3e}",
                    f"sum={np.sum(csv[:,6]):.3e}",
                    f"sum={np.sum(csv[:,7]):.3e}",
                    f"sum={np.sum(csv[:,8]):.3e}",
                    f"sum={np.sum(csv[:,9]):.3e}",
                ]
                times = [csv[i][0] / dic["tscale"] for i in range(csv.shape[0])]
                labels[k] += f" ({fol.split('/')[-1]})"
                axis.step(
                    times,
                    [csv[i][k + 1] for i in range(csv.shape[0])],
                    lw=2,
                    color=dic["colors"][nfol],
                    label=labels[k],
                )
                axis.set_title(plot + f", {dic['case']}")
                axis.set_ylabel(ylabel)
                axis.set_xlabel(f"Time [{dic['tlabel']}]")
                axis.legend()
        dic["fig"].savefig(
            f"{dic['where']}/{dic['case']}_performance{kind}.png", bbox_inches="tight"
        )


def sparse_data(dic):
    """
    Make the plots related to the sparse data (e.g., pressure in sensors over time)

    Args:
        dic (dict): Global dictionary

    Returns:
        None

    """
    dic["fig"] = plt.figure(figsize=(25, 40))
    plots = ["sensors", "boxA", "boxB", "boxC", "facie 1"]
    ylabels = ["Presure [Pa]", "Mass [kg]", "Mass [kg]", "Length [m]", "Mass [kg]"]
    labels = [
        ["p1", "p2"],
        ["mobA", "immA", "dissA", "sealA"],
        ["mobB", "immB", "dissB", "sealB"],
        ["MC"],
        ["sealTot"],
    ]
    dic["nfigs"] = 5
    if dic["case"] != "spe11a":
        plots += ["boundaries"]
        ylabels += ["Mass [kg]"]
        labels.append(["boundTot"])
        dic["nfigs"] += 1
    for k, (plot, ylabel) in enumerate(zip(plots, ylabels)):
        axis = dic["fig"].add_subplot(dic["nfigs"], 3, k + 1)
        for nfol, fol in enumerate(dic["folders"]):
            ncol = sum(len(labels[i]) for i in range(k)) + 1
            axis.text(
                0.7,
                0.15 + nfol * 0.05,
                dic["folders"][-1 - nfol].split("/")[-1],
                transform=axis.transAxes,
                verticalalignment="top",
                bbox=dic["props"],
                color=dic["colors"][len(dic["folders"]) - nfol - 1],
            )
            if os.path.isfile(f"{fol}/data/{dic['case']}_time_series.csv"):
                csv = np.genfromtxt(
                    f"{fol}/data/{dic['case']}_time_series.csv",
                    delimiter=",",
                    skip_header=1,
                )
            else:
                csv = np.genfromtxt(
                    f"{fol}/{dic['case']}_time_series.csv",
                    delimiter=",",
                    skip_header=1,
                )
            times = [csv[i][0] / dic["tscale"] for i in range(csv.shape[0])]
            for j, label in enumerate(labels[k]):
                if nfol == 0:
                    axis.plot(
                        times,
                        [csv[i][ncol] for i in range(csv.shape[0])],
                        label=label,
                        color="k",
                        linestyle=dic["linestyle"][-1 + j],
                    )
                if label == "MC":
                    axis.plot(
                        times,
                        [csv[i][ncol] for i in range(csv.shape[0])],
                        color=dic["colors"][nfol],
                        linestyle=dic["linestyle"][-1 + j + nfol],
                    )
                else:
                    axis.plot(
                        times,
                        [csv[i][ncol] for i in range(csv.shape[0])],
                        color=dic["colors"][nfol],
                        linestyle=dic["linestyle"][-1 + j],
                    )
                ncol += 1
        axis.set_title(plot + f", {dic['case']}")
        axis.set_ylabel(ylabel)
        axis.set_xlabel(f"Time [{dic['tlabel']}]")
        axis.legend()
    dic["fig"].savefig(
        f"{dic['where']}/{dic['case']}_sparse_data.png", bbox_inches="tight"
    )


def generate_grid(dic):
    """
    Create the plotting grid and load the times

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    dic["files"] = [
        f
        for f in os.listdir(f"{dic['folders'][0]}{dic['dataf']}")
        if f.endswith(f"{dic['tlabel']}.csv")
    ]
    dic["times"] = np.array(
        [int(file[19:-5]) for file in dic["files"] if len(file) < 30]
    )
    if dic["times"].size == 0:
        dic["times"] = np.array([int(file[31:-5]) for file in dic["files"]])
    dic["sort_ind"] = np.argsort(dic["times"])
    dic["times"] = [dic["times"][i] for i in dic["sort_ind"]]
    if dic["time"].size == 1:
        if dic["time"] > 0:
            dic["times"] = list(range(0, dic["times"][-1] + 1, dic["time"]))
        else:
            dic["times"] = [dic["time"]]
    else:
        dic["times"] = list(dic["time"])
    csv = np.genfromtxt(
        f"{dic['folders'][0]}{dic['dataf']}/{dic['files'][0]}",
        delimiter=",",
        skip_header=1,
    )
    dic["length"] = csv[-1][0] + csv[0][0]
    dic["width"] = csv[-1][dic["dims"] - 2] + csv[0][dic["dims"] - 2]
    dic["height"] = csv[-1][dic["dims"] - 1] + csv[0][dic["dims"] - 1]
    dic["xmx"] = np.linspace(
        0, dic["length"], round(dic["length"] / (2.0 * csv[0][0])) + 1
    )
    dic["ymy"] = np.linspace(
        0, dic["width"], round(dic["width"] / (2.0 * csv[0][dic["dims"] - 2])) + 1
    )
    dic["zmz"] = np.linspace(
        0, dic["height"], round(dic["height"] / (2.0 * csv[0][dic["dims"] - 1])) + 1
    )
    dic["xmsh"], dic["zmsh"] = np.meshgrid(dic["xmx"], dic["zmz"][::-1])
    if dic["generate"] in ["all", "dense_performance-spatial"]:
        dic["kinds"] = ["", "_performance"]
    elif dic["generate"][:5] == "dense":
        dic["kinds"] = [""]
    else:
        dic["kinds"] = ["_performance"]
    dic["cmaps"] = [
        "seismic",
        "jet",
        "viridis",
        "viridis_r",
        "PuOr",
        "PuOr_r",
        "turbo",
        "coolwarm",
    ]


def handle_kind(dic, kind):
    """
    Identify between dense and performance-spatial

    Args:
        dic (dict): Global dictionary\n
        kind (list): Strings with the type of data to generate

    Returns:
        dic (dict): Modified global dictionary

    """
    if kind == "":
        dic["quantities"] = [
            "pressure",
            "sgas",
            "xco2",
            "xh20",
            "gden",
            "wden",
            "tco2",
        ]
        dic["units"] = [
            "[Pa]",
            "[-]",
            "[-]",
            "[-]",
            r"[kg/m$^3$]",
            r"[kg/m$^3$]",
            "[kg]",
        ]
        dic["allplots"] = [-1, -1, -1, -1, -1, -1, -1]
        if dic["case"] != "spe11a":
            dic["quantities"] += ["temp"]
            dic["units"] += ["C"]
            dic["allplots"] += [-1]
    else:
        dic["quantities"] = [
            "cvol",
            "arat",
            "CO2 max_norm_res",
            "H2O max_norm_res",
            "CO2 mb_error",
            "H2O mb_error",
        ]
        dic["units"] = [r"[m$^3$]", "[-]", "[-]", "[-]", "[-]", "[-]"]
        dic["allplots"] = [0, 0, -1, -1, -1, -1]


def ini_quantity_plot(dic):
    """
    Initialized the size of the Figure according to the spe case

    Args:
        dic (dict): Global dictionary

    Returns:
        dic (dict): Modified global dictionary

    """
    if dic["case"] != "spe11a":
        dic["fig"] = plt.figure(figsize=(50, 3 * len(dic["ptimes"])))
        if dic["lower"]:
            dic["fig"] = plt.figure(figsize=(100, 3 * len(dic["ptimes"])))
    else:
        dic["fig"] = plt.figure(figsize=(45, 6.5 * len(dic["ptimes"])), dpi=80)
    for name in ["plot", "min", "max", "sum"]:
        dic[f"{name}"] = []


def dense_data(dic):
    """
    Make the plots related to the dense data (e.g., saturation maps)

    Args:
        dic (dict): Global dictionary

    Returns:
        None

    """
    generate_grid(dic)
    for kind in dic["kinds"]:
        handle_kind(dic, kind)
        csvs = []
        for tmap in dic["times"]:
            csvs.append(
                np.genfromtxt(
                    f"{dic['folders'][0]}{dic['dataf']}/{dic['case']}{kind}_spatial_map_"
                    + f"{tmap}{dic['tlabel']}.csv",
                    delimiter=",",
                    skip_header=1,
                )
            )
        for k, quantity in enumerate(dic["quantities"]):
            quan = np.array(
                [csvs[0][i][dic["dims"] + k] for i in range(csvs[0].shape[0])]
            )
            if np.isnan(quan).all():
                continue
            dic["minc"], dic["maxc"] = (
                np.min(quan[~np.isnan(quan)]),
                np.max(quan[~np.isnan(quan)]),
            )
            dic["ptimes"] = dic["times"][: dic["allplots"][k]] + [dic["times"][-1]]
            ini_quantity_plot(dic)
            for n, tmap in enumerate(dic["ptimes"]):
                quan = np.array(
                    [csvs[n][i][dic["dims"] + k] for i in range(csvs[n].shape[0])]
                )
                dic["min"].append(np.min(quan[~np.isnan(quan)]))
                dic["max"].append(np.max(quan[~np.isnan(quan)]))
                if quantity == "tco2":
                    dic["sum"].append(np.sum(quan[quan >= 0]))
                dic["minc"] = min(dic["minc"], dic["min"][-1])
                dic["maxc"] = max(dic["maxc"], dic["max"][-1])
                dic["plot"].append(np.zeros([len(dic["zmz"]) - 1, len(dic["xmx"]) - 1]))
                for i in np.arange(0, len(dic["zmz"]) - 1):
                    if dic["case"] != "spe11c":
                        dic["plot"][-1][-1 - i, :] = quan[
                            i * (len(dic["xmx"]) - 1) : (i + 1) * (len(dic["xmx"]) - 1)
                        ]
                    else:
                        dic["plot"][-1][-1 - i, :] = quan[
                            (i * (len(dic["ymy"]) - 1) * (len(dic["xmx"]) - 1))
                            + mt.floor((len(dic["ymy"]) - 1) / 2)
                            * (len(dic["xmx"]) - 1) : (
                                (len(dic["xmx"]) - 1)
                                + i * (len(dic["ymy"]) - 1) * (len(dic["xmx"]) - 1)
                                + mt.floor((len(dic["ymy"]) - 1) / 2)
                                * (len(dic["xmx"]) - 1)
                            )
                        ]
            for j, time in enumerate(dic["ptimes"]):
                if j + 1 < len(dic["ptimes"]):
                    print(
                        f"Plotting {quantity}, time {j+1} out of {len(dic['ptimes'])}",
                        end="\r",
                    )
                else:
                    print(
                        f"Plotting {quantity}, time {j+1} out of {len(dic['ptimes'])}"
                    )
                axis = dic["fig"].add_subplot(len(dic["ptimes"]), 3, j + 1)
                imag = axis.pcolormesh(
                    dic["xmsh"],
                    dic["zmsh"],
                    dic["plot"][j],
                    shading="flat",
                    cmap=dic["cmaps"][k],
                )
                if quantity == "tco2":
                    axis.set_title(
                        f"{time}{dic['tlabel']}, {quantity} "
                        + dic["units"][k]
                        + f"(sum={dic['sum'][j]:.1E})"
                        + f", {dic['case']} ({dic['folders'][0].split('/')[-1]})"
                    )
                else:
                    if dic["allplots"][k] == -1:
                        timet = f"{time}{dic['tlabel']}, "
                    else:
                        timet = ""
                    axis.set_title(
                        timet
                        + f"{quantity} "
                        + dic["units"][k]
                        + f"(min={dic['min'][j]:.1E}, max={dic['max'][j]:.1E})"
                        + f", {dic['case']} ({dic['folders'][0].split('/')[-1]})"
                    )
                axis.axis("scaled")
                axis.xaxis.set_major_locator(ticker.MaxNLocator(14))
                axis.yaxis.set_major_locator(ticker.MaxNLocator(4))
                imag.set_clim(
                    dic["minc"],
                    dic["maxc"],
                )
                if j % 3 != 0:
                    axis.set_yticks([])
                if (
                    j
                    < ((len(dic["ptimes"]) - len(dic["ptimes"]) % 3) / 3) * 3
                    - (3 - len(dic["ptimes"]) % 3)
                    or (len(dic["ptimes"]) % 3 == 1 and j == len(dic["ptimes"]) - 4)
                    or (len(dic["ptimes"]) % 3 == 2 and j == len(dic["ptimes"]) - 5)
                ):
                    axis.set_xticks([])
                if (
                    (j + 1) % 3 == 0
                    or len(dic["ptimes"]) == 1
                    or (len(dic["ptimes"]) == 2 and j == 1)
                ):
                    divider = make_axes_locatable(axis)
                    vect = np.linspace(
                        dic["minc"],
                        dic["maxc"],
                        5,
                        endpoint=True,
                    )
                    dic["fig"].colorbar(
                        imag,
                        cax=divider.append_axes("right", size="5%", pad=0.05),
                        orientation="vertical",
                        ticks=vect,
                        format=lambda x, _: f"{x:.2e}",
                    )
                if dic["lower"]:
                    if dic["case"] == "spe11a":
                        axis.set_ylim([0, 0.55])
                    else:
                        axis.set_ylim([0, 550])
            dic["fig"].savefig(
                f"{dic['where']}/{dic['case']}_{quantity}_2Dmaps.png",
                bbox_inches="tight",
            )
            plt.close()


if __name__ == "__main__":
    main()
