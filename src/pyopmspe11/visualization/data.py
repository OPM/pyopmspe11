# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT
# pylint: disable=C0302, R0912, R0914, R0801, R0915, E1102, C0325

"""
Script to write the benchmark data
"""

import argparse
import warnings
import csv
from io import StringIO
from shapely.geometry import Polygon
from alive_progress import alive_bar
from rtree import index
import numpy as np
from scipy.interpolate import interp1d
from opm.io.ecl import EclFile as OpmFile
from opm.io.ecl import EGrid as OpmGrid
from opm.io.ecl import ERst as OpmRestart
from opm.io.ecl import ESmry as OpmSummary

GAS_DEN_REF = 1.86843
WAT_DEN_REF = 998.108
SECONDS_IN_YEAR = 31536000
SGAS_THR = 0.097


def main():
    """Postprocessing to generate the benchmark data"""
    parser = argparse.ArgumentParser(description="Main script to process the data")
    parser.add_argument(
        "-p",
        "--path",
        default="output",
        help="The name of the output folder ('output' by default).",
    )
    parser.add_argument(
        "-d",
        "--deck",
        default="spe11b",
        help="The simulated case ('spe11b' by default).",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        default="10,1,5",
        help="Number of x, y, and z elements to write the data ('10,1,5' by default).",
    )
    parser.add_argument(
        "-t",
        "--time",
        default="24",
        help="If one number, time step for the spatial maps (spe11a [h]; spe11b/c "
        "[y]) ('24' by default); otherwise, times separated by commas.",
    )
    parser.add_argument(
        "-w",
        "--write",
        default="0.1",
        help="Time interval for the sparse and performance data (spe11a [h]; spe11b/c [y])"
        " ('0.1' by default).",
    )
    parser.add_argument(
        "-g",
        "--generate",
        default="sparse",
        help="Write only the 'dense', 'sparse', 'performance', 'performance-spatial', "
        "'dense_performance', 'performance_sparse', 'dense_performance-spatial', "
        "'dense_sparse', 'dense_performance_sparse', or 'all'",
    )
    parser.add_argument(
        "-n",
        "--neighbourhood",
        default="",
        help="Region to model; valid options are 'lower' or '' (all reservoir) ('' by default)",
    )
    parser.add_argument(
        "-f",
        "--subfolders",
        default=1,
        help="Set to 0 to not create the subfolders deck, flow, data, and figures, i.e., to "
        "write all generated files in the output directory ('1' by default).",
    )
    parser.add_argument(
        "-s",
        "--showpywarn",
        default=0,
        help="Set to 1 to show Python warnings ('0' by default).",
    )
    cmdargs = vars(parser.parse_known_args()[0])
    if int(cmdargs["showpywarn"]) != 1:  # Show or hidde python warnings
        warnings.warn = lambda *args, **kwargs: None
    dig = {"path": cmdargs["path"].strip()}
    dig["case"] = cmdargs["deck"].strip()
    dig["mode"] = cmdargs["generate"].strip()
    dig["lower"] = bool(cmdargs["neighbourhood"].strip())  # Lower model
    if int(cmdargs["subfolders"]) == 1:
        dig["deckf"] = f"{dig['path']}/deck"
        dig["flowf"] = f"{dig['path']}/flow"
        dig["where"] = f"{dig['path']}/data"
    else:
        dig["deckf"] = dig["path"]
        dig["flowf"] = dig["path"]
        dig["where"] = dig["path"]
    dig["nxyz"] = np.genfromtxt(
        StringIO(cmdargs["resolution"]), delimiter=",", dtype=int
    )
    dig["sim"] = dig["flowf"] + f"/{dig['path'].split('/')[-1].upper()}"
    if dig["case"] == "spe11a":
        dig["dense_t"] = (
            np.genfromtxt(StringIO(cmdargs["time"]), delimiter=",", dtype=float) * 3600
        )
        dig["sparse_t"] = 1.0 * round(float(cmdargs["write"].strip()) * 3600)
        dig["dims"] = [2.8, 1.0, 1.2]
        dig["dof"], dig["nxyz"][1] = 2, 1
    else:
        dig["dense_t"] = (
            np.genfromtxt(StringIO(cmdargs["time"]), delimiter=",", dtype=float)
            * SECONDS_IN_YEAR
        )
        dig["sparse_t"] = float(cmdargs["write"].strip()) * SECONDS_IN_YEAR
        dig["dims"] = [8400.0, 1.0, 1200.0]
        dig["dof"] = 3
    if dig["case"] == "spe11c":
        dig["dims"][1] = 5000.0
    dig["nocellsr"] = dig["nxyz"][0] * dig["nxyz"][1] * dig["nxyz"][2]
    dig["noxzr"] = dig["nxyz"][0] * dig["nxyz"][2]
    dig["time_initial"], dig["no_skip_rst"], dig["times"] = 0, 0, []
    dig["immiscible"] = False
    read_simulations(dig)
    if dig["mode"] in [
        "performance",
        "all",
        "dense_performance",
        "performance_sparse",
        "dense_performance_sparse",
    ]:
        performance(dig)
    if dig["mode"] in [
        "all",
        "sparse",
        "dense_sparse",
        "dense_performance_sparse",
        "performance_sparse",
    ]:
        sparse_data(dig)
    if dig["mode"] in [
        "all",
        "performance-spatial",
        "dense",
        "dense_performance",
        "dense_sparse",
        "dense_performance-spatial",
        "dense_performance_sparse",
    ]:
        if isinstance(dig["dense_t"], float):
            dig["dense_t"] = [
                i * dig["dense_t"]
                for i in range(int(np.floor((dig["times"][-1]) / dig["dense_t"])) + 1)
            ]
        dense_data(dig)
    print(f"The csv files have been written to {dig['where']}")


def read_simulations(dig):
    """
    Read the simulation files using OPM

    Args:
        dig (dict): Global dictionary

    Returns:
        dig (dict): Modified global dictionary

    """
    dig["unrst"] = OpmRestart(f"{dig['sim']}.UNRST")
    time = []
    for i in range(len(dig["unrst"].report_steps)):
        time.append(86400 * dig["unrst"]["DOUBHEAD", i][0])
        if len(dig["times"]) == 0:
            if dig["unrst"].count("RSW", 0):
                if np.max(dig["unrst"]["RSW", i]) > 0:
                    dig["time_initial"] = 86400 * dig["unrst"]["DOUBHEAD", i - 1][0]
                    dig["no_skip_rst"] = i - 1
                    dig["times"].append(0)
                    dig["times"].append(time[-1] - dig["time_initial"])
            else:
                dig["immiscible"] = True
                if np.max(dig["unrst"]["SGAS", i]) > 0:
                    dig["time_initial"] = 86400 * dig["unrst"]["DOUBHEAD", i - 1][0]
                    dig["no_skip_rst"] = i - 1
                    dig["times"].append(0)
                    dig["times"].append(time[-1] - dig["time_initial"])
        else:
            dig["times"].append(time[-1] - dig["time_initial"])
    if not dig["times"]:
        dig["times"] = time
    dig["init"] = OpmFile(f"{dig['sim']}.INIT")
    dig["egrid"] = OpmGrid(f"{dig['sim']}.EGRID")
    dig["smspec"] = OpmSummary(f"{dig['sim']}.SMSPEC")
    dig["norst"] = len(dig["unrst"].report_steps)
    dig["porv"] = np.array(dig["init"]["PORV"])
    dig["actind"] = list(i for i, porv in enumerate(dig["porv"]) if porv > 0)
    dig["porva"] = np.array([porv for porv in dig["porv"] if porv > 0])
    dig["nocellst"], dig["nocellsa"] = (
        len(dig["porv"]),
        dig["egrid"].active_cells,
    )
    dig["times_summary"] = [0.0]
    dig["times_summary"] += list(86400.0 * dig["smspec"]["TIME"] - dig["time_initial"])
    dig["gxyz"] = [
        dig["egrid"].dimension[0],
        dig["egrid"].dimension[1],
        dig["egrid"].dimension[2],
    ]
    dig["noxz"] = dig["egrid"].dimension[0] * dig["egrid"].dimension[2]
    dig["cp"] = dig["porv"][-1] == 0


def performance(dig):
    """
    Generate the performance within the benchmark format

    Args:
        dig (dict): Global dictionary

    Returns:
        None

    """
    dil = {"infosteps": []}
    dil["times_data"] = np.linspace(
        0, dig["times"][-1], round(dig["times"][-1] / dig["sparse_t"]) + 1
    )
    with open(
        f"{dig['flowf']}/{dig['path'].split('/')[-1].upper()}.INFOSTEP",
        "r",
        encoding="utf8",
    ) as file:
        for j, row in enumerate(csv.reader(file)):
            if j == 0:
                tag = (row[0].strip()).split()
            else:
                if float((row[0].strip()).split()[0]) >= (
                    (dig["time_initial"] - dig["sparse_t"]) / 86400.0
                ):
                    dil["infosteps"].append(
                        [float(column) for column in (row[0].strip()).split()]
                    )
    infotimes = [
        infostep[tag.index("Time(day)")] * 86400 - dig["time_initial"]
        for infostep in dil["infosteps"]
    ]
    time0 = max(0, dig["no_skip_rst"] - 1)
    dil["map_info"] = np.array(
        [time0 + int(np.floor(infotime / dig["sparse_t"])) for infotime in infotimes]
    )
    dil["detail_info"] = [0]
    for j, infotime in enumerate(infotimes):
        if j < len(infotimes) - 1:
            if infotime == infotimes[j + 1]:
                dil["detail_info"].append(dil["detail_info"][-1])
            else:
                dil["detail_info"].append(1 + dil["detail_info"][-1])
    dil["detail_info"] = np.array(dil["detail_info"])
    dil["times_det"] = []
    infotimes = np.array(infotimes)
    for j in range(np.max(dil["detail_info"]) + 1):
        ind = j == dil["detail_info"]
        dil["times_det"].append(np.max(infotimes[ind]))
    dil["times_det"] = np.array(dil["times_det"])
    dil["fsteps"] = np.array(
        [1.0 * (infostep[tag.index("Conv")] == 0) for infostep in dil["infosteps"]]
    )
    dil["nress"] = np.array(
        [infostep[tag.index("Lins")] for infostep in dil["infosteps"]]
    )
    dil["tlinsols"] = np.array(
        [infostep[tag.index("LSolve")] for infostep in dil["infosteps"]]
    )
    dil["liniters"] = np.array(
        [infostep[tag.index("LinIt")] for infostep in dil["infosteps"]]
    )
    dil["nliters"] = np.array(
        [infostep[tag.index("NewtIt")] for infostep in dil["infosteps"]]
    )
    dil["tsteps"] = np.array(
        [
            86400 * infostep[tag.index("TStep(day)")] * infostep[tag.index("Conv")]
            for infostep in dil["infosteps"]
        ]
    )
    dil["alltsteps"] = np.array(
        [86400 * infostep[tag.index("TStep(day)")] for infostep in dil["infosteps"]]
    )
    tcpu = dig["smspec"]["TCPU"]
    fgmip = dig["smspec"]["FGMIP"]
    times = 86400.0 * dig["smspec"]["TIME"] - dig["time_initial"]
    dil["map_sum"] = np.array(
        [time0 + int(np.floor(time / dig["sparse_t"])) for time in dil["times_det"]]
    )
    if dig["time_initial"] > 0:
        tcpu = tcpu[-len(dil["map_sum"]) - 1 :]
        tcpu = tcpu[1:] - tcpu[:-1]
    else:
        tcpu = tcpu[-len(dil["map_sum"]) :]
        tcpu[1:] -= tcpu[:-1]
    if dig["time_initial"] == 0:
        times = np.insert(times, 0, 0)
        fgmip = np.insert(fgmip, 0, 0)
    interp_fgmip = interp1d(
        times,
        fgmip,
        fill_value="extrapolate",
    )
    write_performance(dig, dil, interp_fgmip, tcpu, infotimes)


def write_performance(dig, dil, interp_fgmip, tcpu, infotimes):
    """
    Write the performance data

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary\n
        interp_fgmip (object): Interpolator (time) for the CO2 mass\n
        tcpu (array): Floats with the simulation times\n
        infotimes (array): Floats with the simulation time steps

    Returns:
        None

    """
    dil["text"] = []
    dil["text"].append(
        "# t [s], tstep [s], fsteps [-], mass [kg], dof [-], nliter [-], "
        + "nres [-], liniter [-], runtime [s], tlinsol [s]"
    )
    if dig["no_skip_rst"] == 0:
        dil["text"].append(
            "0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, "
            + f"{dig['dof'] * dig['nocellsa']:.3e}, "
            + "0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00"
        )
        dil["times_data"] = np.delete(dil["times_data"], 0)
    freq = [0]
    for j, time in enumerate(dil["times_data"]):
        itd = j == dil["map_sum"]
        if np.sum(tcpu[itd]) == 0:
            freq.append(freq[-1] + 1)
        else:
            freq.append(0)
    freq.pop(0)
    weig = []
    quan = 1
    for val in freq[::-1]:
        if val > 0 and quan == 1:
            quan = val + 1.0
        elif val == 0:
            weig.append(quan)
            quan = 1.0
            continue
        weig.append(quan)
    weig = weig[::-1]
    for j, time in enumerate(dil["times_data"]):
        if freq[j] == 0:
            ind = j == dil["map_info"]
            itd = j == dil["map_sum"]
            dil["tstep"] = np.sum(dil["tsteps"][ind])
            if dil["tstep"] > 0:
                dil["tstep"] /= np.sum(ind)
        dil["text"].append(
            f"{time:.3e}, "
            + f"{dil['tstep']/weig[j]:.3e}, "
            + f"{np.sum(dil['fsteps'][ind])/weig[j]:.3e}, "
            + f"{interp_fgmip(time):.3e}, "
            + f"{dig['dof'] * dig['nocellsa']:.3e}, "
            + f"{np.sum(dil['nliters'][ind])/weig[j]:.3e}, "
            + f"{np.sum(dil['nress'][ind])/weig[j]:.3e}, "
            + f"{np.sum(dil['liniters'][ind])/weig[j]:.3e}, "
            + f"{np.sum(tcpu[itd])/weig[j]:.3e}, "
            + f"{np.sum(dil['tlinsols'][ind])/weig[j]:.3e}"
        )
    with open(
        f"{dig['where']}/{dig['case']}_performance_time_series.csv",
        "w",
        encoding="utf8",
    ) as file:
        file.write("\n".join(dil["text"]))
    dil["text"] = []
    dil["text"].append(
        "# t [s], tstep [s], fsteps [-], mass [kg], dof [-], nliter [-], "
        + "nres [-], liniter [-], runtime [s], tlinsol [s]"
    )
    for j in range(np.max(dil["detail_info"]) + 1):
        ind = j == dil["detail_info"]
        time = np.max(infotimes[ind])
        if time >= 0:
            dil["text"].append(
                f"{time:.3e}, "
                + f"{np.max(dil['tsteps'][ind]):.3e}, "
                + f"{np.sum(dil['fsteps'][ind]):.3e}, "
                + f"{interp_fgmip(time):.3e}, "
                + f"{dig['dof'] * dig['nocellsa']:.3e}, "
                + f"{np.sum(dil['nliters'][ind]):.3e}, "
                + f"{np.sum(dil['nress'][ind]):.3e}, "
                + f"{np.sum(dil['liniters'][ind]):.3e}, "
                + f"{tcpu[j]:.3e}, "
                + f"{np.sum(dil['tlinsols'][ind]):.3e}"
            )
    with open(
        f"{dig['where']}/{dig['case']}_performance_time_series_detailed.csv",
        "w",
        encoding="utf8",
    ) as file:
        file.write("\n".join(dil["text"]))


def create_from_summary(dig, dil):
    """
    Use the summary arrays for the sparse data interpolation

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary

    Returns:
        dil (dict): Modified local dictionary

    """
    ind, names, i_jk = 0, [], []
    for key in dig["smspec"].keys():
        if key[0 : len("BWPR")] == "BWPR" and "," in key[len("BWPR") + 1 :]:
            names.append(key)
            i_jk.append(
                np.genfromtxt(
                    StringIO(key[len("BWPR") + 1 :]), delimiter=",", dtype=int
                )[0]
            )
            ind += 1
            if ind == 2:
                break
    sort = sorted(range(len(i_jk)), key=i_jk.__getitem__)
    pop1 = (
        dig["unrst"]["PRESSURE", 0][dil["fipnum"].index(8)]
        - dig["unrst"]["PCGW", 0][dil["fipnum"].index(8)]
    )
    dil["pop1"] = [pop1 * 1.0e5] + list(dig["smspec"][names[sort[0]]] * 1.0e5)  # Pa
    if not dig["lower"]:
        pop2 = (
            dig["unrst"]["PRESSURE", 0][dil["fipnum"].index(9)]
            - dig["unrst"]["PCGW", 0][dil["fipnum"].index(9)]
        )
        dil["pop2"] = [pop2 * 1.0e5] + list(dig["smspec"][names[sort[1]]] * 1.0e5)  # Pa
    else:
        dil["pop2"] = dil["pop1"]
    for i in dil["fip_diss_a"]:
        dil["moba"] += dig["smspec"][f"RGKMO:{i}"]
        dil["imma"] += dig["smspec"][f"RGKTR:{i}"]
        dil["dissa"] += dig["smspec"][f"RGMDS:{i}"]
    for i in dil["fip_seal_a"]:
        dil["seala"] += (
            dig["smspec"][f"RGMDS:{i}"]
            + dig["smspec"][f"RGKMO:{i}"]
            + dig["smspec"][f"RGKTR:{i}"]
        )
    for i in dil["fip_diss_b"]:
        dil["mobb"] += dig["smspec"][f"RGKMO:{i}"]
        dil["immb"] += dig["smspec"][f"RGKTR:{i}"]
        dil["dissb"] += dig["smspec"][f"RGMDS:{i}"]
    for i in dil["fip_seal_b"]:
        dil["sealb"] += (
            dig["smspec"][f"RGMDS:{i}"]
            + dig["smspec"][f"RGKMO:{i}"]
            + dig["smspec"][f"RGKTR:{i}"]
        )
    dil["sealt"] = dil["seala"] + dil["sealb"]
    for name in ["RGMDS", "RGKMO", "RGKTR"]:
        if not dig["lower"]:
            dil["sealt"] += dig["smspec"][f"{name}:7"] + dig["smspec"][f"{name}:9"]
        else:
            if 7 in dil["fipnum"]:
                dil["sealt"] += dig["smspec"][f"{name}:7"]
    if dig["case"] != "spe11a":
        if 10 in dil["fipnum"]:
            sealbound = (
                dig["smspec"]["RGMDS:10"]
                + dig["smspec"]["RGKMO:10"]
                + dig["smspec"]["RGKTR:10"]
            )
            dil["sealt"] += sealbound
            dil["boundtot"] = sealbound
        else:
            dil["boundtot"] = 0
        for i in dil["fip_bound_t"]:
            if i in dil["fipnum"]:
                dil["boundtot"] += (
                    dig["smspec"][f"RGMDS:{i}"]
                    + dig["smspec"][f"RGKMO:{i}"]
                    + dig["smspec"][f"RGKTR:{i}"]
                )


def sparse_data(dig):
    """
    Generate the sparse data within the benchmark format

    Args:
        dig (dict): Global dictionary

    Returns:
        None

    """
    dil = {
        "times_data": np.linspace(
            0, dig["times"][-1], round(dig["times"][-1] / dig["sparse_t"]) + 1
        )
    }
    dil["fipnum"] = list(dig["init"]["FIPNUM"])
    for name in ["dx", "dy", "dz"]:
        dil[f"{name}"] = np.array(dig["init"][name.upper()])
    dil["names"] = [
        "pop1",
        "pop2",
        "moba",
        "imma",
        "dissa",
        "seala",
        "mobb",
        "immb",
        "dissb",
        "sealb",
        "sealt",
    ]
    if dig["case"] != "spe11a":
        dil["names"] += ["boundtot"]
    for ent in dil["names"]:
        dil[ent] = 0.0
    dil["m_c"] = []
    handle_fipnums(dig, dil)
    create_from_summary(dig, dil)
    # Using the restart data
    if not dig["immiscible"]:
        compute_m_c(dig, dil)
    else:
        dil["m_c"] = [0.0] * (dig["norst"] - dig["no_skip_rst"] - 1)
    if dig["lower"]:
        for name in ["seala", "mobb", "immb", "dissb", "sealb"]:
            dil[name] = [0.0] * (len(dig["times_summary"]) - 1)
    write_sparse_data(dig, dil)


def handle_fipnums(dig, dil):
    """
    Set the fipnum groups to compute the sparse data

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary

    Returns:
        dil (dict): Modified local dictionary

    """
    if dig["lower"]:
        dil["fip_diss_a"] = [2, 4, 8]
        dil["fip_seal_a"] = []
        dil["fip_diss_b"] = []
        dil["fip_seal_b"] = []
    else:
        dil["fip_diss_a"] = [2, 4, 5, 8]
        dil["fip_seal_a"] = [5, 8]
        dil["fip_diss_b"] = [3, 6]
        dil["fip_seal_b"] = [6]
    if dig["case"] != "spe11a":
        dil["fip_bound_t"] = [11]
    if dig["case"] == "spe11c":
        dil["fip_diss_a"] += [13, 17]
        dil["fip_bound_t"] += [13, 17]
        if not dig["lower"]:
            dil["fip_diss_a"] += [14]
            dil["fip_seal_a"] += [14]
            dil["fip_diss_b"] += [15, 16]
            dil["fip_seal_b"] += [16]
            dil["fip_bound_t"] += [14, 15, 16]
        if 18 in dil["fipnum"]:
            dil["fip_diss_a"] += [18]
            dil["fip_seal_a"] += [18]
            dil["fip_bound_t"] += [18]
        if 12 in dil["fipnum"]:
            dil["fip_diss_a"] += [12]
            dil["fip_seal_a"] += [12]


def compute_m_c(dig, dil):
    """
    Normalized total variation of the concentration field within Box C

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary

    Returns:
        dil (dict): Modified local dictionary

    """
    dil["boxc"] = np.array([fip in (4, 12, 17, 18) for fip in dil["fipnum"]])
    dil["boxc_x"] = np.roll(dil["boxc"], 1)
    dil["boxc_y"] = np.roll(dil["boxc"], -dig["gxyz"][0])
    dil["boxc_z"] = np.roll(dil["boxc"], -dig["gxyz"][0] * dig["gxyz"][1])
    for t_n in range(dig["no_skip_rst"] + 1, dig["norst"]):
        rss = np.array(dig["unrst"]["RSW", t_n])
        dil["xcw"] = np.divide(rss, rss + WAT_DEN_REF / GAS_DEN_REF)
        rssat = np.array(dig["unrst"]["RSWSAT", t_n])
        x_l_co2_max = np.divide(rssat, rssat + WAT_DEN_REF / GAS_DEN_REF)
        dil["xcw"] = np.divide(dil["xcw"], x_l_co2_max)
        if dig["case"] != "spe11c":
            dil["m_c"].append(
                np.sum(
                    np.abs(
                        (dil["xcw"][dil["boxc_x"]] - dil["xcw"][dil["boxc"]])
                        * dil["dz"][dil["boxc"]]
                    )
                    + np.abs(
                        (dil["xcw"][dil["boxc_z"]] - dil["xcw"][dil["boxc"]])
                        * dil["dx"][dil["boxc"]]
                    )
                )
            )
        else:
            dil["m_c"].append(
                np.sum(
                    np.abs(
                        (dil["xcw"][dil["boxc_x"]] - dil["xcw"][dil["boxc"]])
                        * dil["dy"][dil["boxc"]]
                        * dil["dz"][dil["boxc"]]
                    )
                    + np.abs(
                        (dil["xcw"][dil["boxc_y"]] - dil["xcw"][dil["boxc"]])
                        * dil["dx"][dil["boxc"]]
                        * dil["dz"][dil["boxc"]]
                    )
                    + np.abs(
                        (dil["xcw"][dil["boxc_z"]] - dil["xcw"][dil["boxc"]])
                        * dil["dx"][dil["boxc"]]
                        * dil["dy"][dil["boxc"]]
                    )
                )
            )


def write_sparse_data(dig, dil):
    """
    Write the sparse data

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary

    Returns:
        None

    """
    for name in dil["names"] + ["m_c"]:
        if name == "m_c":
            interp = interp1d(
                dig["times"],
                [0.0] + dil[f"{name}"],
                fill_value="extrapolate",
            )
        elif "pop" in name:
            interp = interp1d(
                dig["times_summary"],
                dil[f"{name}"],
                fill_value="extrapolate",
            )
        else:
            if isinstance(dil[f"{name}"], float):
                dil[f"{name}"] = [0.0] * (len(dig["times_summary"]) - 1)
            interp = interp1d(
                dig["times_summary"],
                [0.0] + list(dil[f"{name}"]),
                fill_value="extrapolate",
            )
        dil[f"{name}"] = interp(dil["times_data"])
    text = [
        "# t [s], p1 [Pa], p2 [Pa], mobA [kg], immA [kg], dissA [kg], sealA [kg], "
        + "mobB [kg], immB [kg], dissB [kg], sealB [kg], MC [m], sealTot [kg]"
    ]
    if dig["case"] == "spe11a":
        for j, time in enumerate(dil["times_data"]):
            text.append(
                f"{time:.3e}, {dil['pop1'][j]:.5e}, {dil['pop2'][j]:.5e}, "
                f"{dil['moba'][j]:.3e}, {dil['imma'][j]:.3e}, {dil['dissa'][j]:.3e}, "
                f"{dil['seala'][j]:.3e}, {dil['mobb'][j]:.3e}, {dil['immb'][j]:.3e}, "
                f"{dil['dissb'][j]:.3e}, {dil['sealb'][j]:.3e}, {dil['m_c'][j]:.3e}, "
                f"{dil['sealt'][j]:.3e}"
            )
    else:
        text[-1] += ", boundTot [kg]"
        for j, time in enumerate(dil["times_data"]):
            text.append(
                f"{time:.4e}, {dil['pop1'][j]:.3e}, {dil['pop2'][j]:.3e}, "
                f"{dil['moba'][j]:.3e}, {dil['imma'][j]:.3e}, {dil['dissa'][j]:.3e}, "
                f"{dil['seala'][j]:.3e}, {dil['mobb'][j]:.3e}, {dil['immb'][j]:.3e}, "
                f"{dil['dissb'][j]:.3e}, {dil['sealb'][j]:.3e}, {dil['m_c'][j]:.3e}, "
                f"{dil['sealt'][j]:.3e}, {dil['boundtot'][j]:.3e}"
            )
    with open(
        f"{dig['where']}/{dig['case']}_time_series.csv", "w", encoding="utf8"
    ) as file:
        file.write("\n".join(text))


def get_corners(dig, dil):
    """
    Get the cell corners from the simulation grid

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary

    Returns:
        dil (dict): Modified local dictionary

    """
    for i in ["x", "z"]:
        dil[f"sim{i}cent"] = [0.0] * dig["noxz"]
    dil["simycent"] = [0.0] * dig["gxyz"][1]
    dil["simpoly"] = [0.0] * dig["noxz"]
    z_0 = 155.04166666666666 if dig["case"] == "spe11c" else 0
    dil["satnum"] = list(dig["init"]["SATNUM"])
    for j in range(dig["gxyz"][2]):
        for i in range(dig["gxyz"][0]):
            n = i + (dig["gxyz"][2] - j - 1) * dig["gxyz"][0]
            xyz = dig["egrid"].xyz_from_ijk(i, 0, dig["gxyz"][2] - j - 1)
            dil["simpoly"][n] = Polygon(
                [
                    [xyz[0][0], dig["dims"][2] - (xyz[2][0] - z_0)],
                    [xyz[0][1], dig["dims"][2] - (xyz[2][1] - z_0)],
                    [xyz[0][5], dig["dims"][2] - (xyz[2][5] - z_0)],
                    [xyz[0][4], dig["dims"][2] - (xyz[2][4] - z_0)],
                ]
            )
            pxz = dil["simpoly"][n].centroid.wkt
            pxz = list(float(j) for j in pxz[7:-1].split(" "))
            dil["simxcent"][n] = pxz[0]
            dil["simzcent"][n] = pxz[1]
            if (
                pxz[1] / 1.2 < 1e-4 / dig["dims"][2]
                and 2.72 / 2.8 < pxz[0] / dig["dims"][0]
            ):
                dil["simxcent"][n] = -1e10
                dil["simzcent"][n] = -1e10
                dil["simpoly"][n] = Polygon(
                    [
                        [xyz[0][0], dig["dims"][2] - (xyz[2][0] - z_0)],
                        [xyz[0][1], dig["dims"][2] - (xyz[2][1] - z_0)],
                        [xyz[0][1], dig["dims"][2] - (xyz[2][1] - z_0)],
                        [xyz[0][0], dig["dims"][2] - (xyz[2][0] - z_0)],
                    ]
                )
    for j in range(dig["gxyz"][1]):
        xyz = dig["egrid"].xyz_from_ijk(0, j, 0)
        dil["simycent"][j] = 0.5 * (xyz[1][2] - xyz[1][1]) + xyz[1][1]
    dil["simycent"] = np.array(dil["simycent"])
    if dig["lower"] and dig["cp"]:
        dil["simzcent"] = np.array(dil["simzcent"])
        dil["simxcent"] = np.insert(
            dil["simxcent"], 0, dil["simxcent"][: dig["gxyz"][0]]
        )
        dil["simzcent"] = np.insert(
            dil["simzcent"], 0, dil["simzcent"][: dig["gxyz"][0]] + 1e-4
        )


def dense_data(dig):
    """
    Generate the dense data within the benchmark format

    Args:
        dig (dict): Global dictionary

    Returns:
        None

    """
    dil = {"rstno": []}
    for time in dig["dense_t"]:
        dil["rstno"].append(dig["times"].index(time))
    get_corners(dig, dil)
    dil["nrstno"] = len(dil["rstno"])
    for i, j, k in zip(["x", "y", "z"], dig["dims"], dig["nxyz"]):
        dil[f"ref{i}vert"] = np.linspace(0, j, k + 1)
        dil[f"ref{i}cent"] = 0.5 * (dil[f"ref{i}vert"][1:] + dil[f"ref{i}vert"][:-1])
    dil["refxgrid"] = np.zeros(dig["nxyz"][0] * dig["nxyz"][2])
    dil["refzgrid"] = np.zeros(dig["nxyz"][0] * dig["nxyz"][2])
    dil["refpoly"] = []
    dil["cell_ind"] = [[] for _ in range(dig["noxz"])]
    dil["cell_indc"] = np.zeros(dig["noxz"], dtype=int)
    dil["cell_cent"] = np.zeros(dig["noxzr"], dtype=float)
    dx = dig["init"]["DX"][: dig["gxyz"][0]]
    dz = dig["init"]["DZ"]
    iszunif = np.min(dz) == np.max(dz) and not dig["lower"]
    if (
        iszunif
        and dig["nxyz"][2] == dig["gxyz"][2]
        and np.min(dx) == np.max(dx)
        and dig["nxyz"][0] == dig["gxyz"][0]
    ):
        for k in range(dig["nxyz"][2]):
            dil["cell_cent"][
                (dig["nxyz"][2] - k - 1)
                * dig["nxyz"][0] : (dig["nxyz"][2] - k)
                * dig["nxyz"][0]
            ] = range(k * dig["nxyz"][0], (k + 1) * dig["nxyz"][0])
        dil["cell_indc"] = dil["cell_cent"]
        for n, row in enumerate(dil["cell_indc"]):
            dil["cell_ind"][n] = [[int(row), 1]]
    elif (
        iszunif
        and dig["nxyz"][2] == dig["gxyz"][2]
        and np.min(dx[2:-2]) == np.max(dx[2:-2])
        and dig["nxyz"][0] == dig["gxyz"][0] - 2
    ):
        for k in range(dig["nxyz"][2]):
            dil["cell_indc"][(dig["gxyz"][2] - k - 1) * dig["gxyz"][0]] = (
                k * dig["nxyz"][0]
            )
            dil["cell_indc"][
                (dig["gxyz"][2] - k - 1) * dig["gxyz"][0]
                + 1 : (dig["gxyz"][2] - k) * dig["gxyz"][0]
            ] = range(k * dig["nxyz"][0], (k + 1) * dig["nxyz"][0] + 1)
            dil["cell_indc"][(dig["gxyz"][2] - k) * dig["gxyz"][0] - 1] = (k + 1) * dig[
                "nxyz"
            ][0] - 1
            dil["cell_cent"][
                (dig["nxyz"][2] - k - 1)
                * dig["nxyz"][0] : (dig["nxyz"][2] - k)
                * dig["nxyz"][0]
            ] = [
                row + 2 * k
                for row in range(k * dig["nxyz"][0] + 1, (k + 1) * dig["nxyz"][0] + 1)
            ]
        for n, row in enumerate(dil["cell_indc"]):
            dil["cell_ind"][n] = [[int(row), 1]]
    elif (
        iszunif
        and dig["gxyz"][2] % dig["nxyz"][2] == 0
        and np.min(dx) == np.max(dx)
        and dig["gxyz"][0] % dig["nxyz"][0] == 0
    ):
        x_n = int(dig["gxyz"][0] / dig["nxyz"][0])
        z_n = int(dig["gxyz"][2] / dig["nxyz"][2])
        for k in range(dig["nxyz"][2]):
            dil["cell_cent"][
                (dig["nxyz"][2] - k - 1)
                * dig["nxyz"][0] : (dig["nxyz"][2] - k)
                * dig["nxyz"][0]
            ] = [
                row * x_n
                + (x_n / 2 - 1)
                + (z_n / 2 - 1) * dig["gxyz"][0]
                + k * (z_n - 1) * dig["gxyz"][0]
                for row in range(k * dig["nxyz"][0], (k + 1) * dig["nxyz"][0])
            ]
        for k in range(dig["nxyz"][2]):
            for k_n in range(z_n):
                for i in range(dig["nxyz"][0]):
                    dil["cell_indc"][
                        (dig["gxyz"][2] - (k * (z_n) + k_n) - 1) * dig["gxyz"][0]
                        + i
                        * (x_n) : (dig["gxyz"][2] - (k * (z_n) + k_n) - 1)
                        * dig["gxyz"][0]
                        + (i + 1) * (x_n)
                    ] = [i + k * dig["nxyz"][0] for _ in range(0, x_n)]
        for n, row in enumerate(dil["cell_indc"]):
            dil["cell_ind"][n] = [[int(row), 1]]
    elif (
        iszunif
        and dig["gxyz"][2] % dig["nxyz"][2] == 0
        and np.min(dx[2:-2]) == np.max(dx[2:-2])
        and (dig["gxyz"][0] - 2) % dig["nxyz"][0] == 0
    ):
        x_n = int((dig["gxyz"][0] - 2) / dig["nxyz"][0])
        z_n = int(dig["gxyz"][2] / dig["nxyz"][2])
        for k in range(dig["nxyz"][2]):
            dil["cell_cent"][
                (dig["nxyz"][2] - k - 1)
                * dig["nxyz"][0] : (dig["nxyz"][2] - k)
                * dig["nxyz"][0]
            ] = [
                row * x_n
                + (x_n / 2 - 1)
                + (z_n / 2 - 1) * dig["gxyz"][0]
                + k * (z_n - 1) * dig["gxyz"][0]
                + 2 * k
                + 1
                for row in range(k * dig["nxyz"][0], (k + 1) * dig["nxyz"][0])
            ]
        for k in range(dig["nxyz"][2]):
            for k_n in range(z_n):
                for i in range(dig["nxyz"][0]):
                    dil["cell_indc"][
                        (dig["gxyz"][2] - (k * (z_n) + k_n) - 1) * dig["gxyz"][0]
                        + i * (x_n)
                    ] = (k * dig["nxyz"][0] + i)
                    dil["cell_indc"][
                        (dig["gxyz"][2] - (k * (z_n) + k_n) - 1) * dig["gxyz"][0]
                        + i * (x_n)
                        + 1 : (dig["gxyz"][2] - (k * (z_n) + k_n) - 1) * dig["gxyz"][0]
                        + (i + 1) * (x_n)
                        + 1
                    ] = [i + k * dig["nxyz"][0] for _ in range(0, x_n)]
                    dil["cell_indc"][
                        (dig["gxyz"][2] - (k * (z_n) + k_n) - 1) * dig["gxyz"][0]
                        + (i + 1) * (x_n)
                        + 1
                    ] = (i + k * dig["nxyz"][0])
        for n, row in enumerate(dil["cell_indc"]):
            dil["cell_ind"][n] = [[int(row), 1]]
    else:
        handle_find_cells_ids(dil)
    dig["actindr"] = []
    if np.max(dil["satnum"]) < 7 and dig["case"] == "spe11a":
        handle_inactive_mapping(dig, dil)
    if dig["case"] == "spe11c":
        handle_yaxis_mapping_intensive(dig, dil)
        handle_yaxis_mapping_extensive(dig, dil)
    if dig["mode"] == "all" or dig["mode"][:5] == "dense":
        names = ["pressure", "sgas", "xco2", "xh20", "gden", "wden", "tco2"]
        if dig["case"] != "spe11a":
            names = ["temp"] + names
        for i, rst in enumerate(dil["rstno"]):
            if i + 1 < dil["nrstno"]:
                print(f"Processing dense data {i+1} out of {dil['nrstno']}", end="\r")
            else:
                print(f"Processing dense data {i+1} out of {dil['nrstno']}")
            t_n = rst + dig["no_skip_rst"]
            generate_arrays(dig, dil, names, t_n)
            map_to_report_grid(dig, dil, names)
            write_dense_data(dig, dil, i)
    if dig["mode"] in ["all", "performance-spatial", "dense_performance-spatial"]:
        handle_performance_spatial(dig, dil)


def handle_find_cells_ids(dil):
    """
    Find the cells ids when reporting and simulation grids differ

    Args:
        dil (dict): Local dictionary

    Returns:
        dil (dict): Modified local dictionary

    """
    ind, idx = 0, index.Index()
    for k, zcen in enumerate(dil["refzcent"]):
        for i, xcen in enumerate(dil["refxcent"]):
            dil["refxgrid"][ind] = xcen
            dil["refzgrid"][ind] = zcen
            dil["refpoly"].append(
                Polygon(
                    [
                        (dil["refxvert"][i], dil["refzvert"][k]),
                        (dil["refxvert"][i + 1], dil["refzvert"][k]),
                        (dil["refxvert"][i + 1], dil["refzvert"][k + 1]),
                        (dil["refxvert"][i], dil["refzvert"][k + 1]),
                    ]
                )
            )
            idx.insert(ind, dil["refpoly"][-1].bounds)
            ind += 1
    print("Processing polygon intersections between simulation and reporting grids")
    with alive_bar(len(dil["simpoly"])) as bar_animation:
        for k, simp in enumerate(dil["simpoly"]):
            bar_animation()
            ovrl = list(idx.intersection(simp.bounds))
            if simp.area > 0:
                for ind in ovrl:
                    area = simp.intersection(dil["refpoly"][ind]).area / simp.area
                    if area > 0:
                        dil["cell_ind"][k].append([ind, area])
                        dil["cell_indc"][k] = ind
            else:
                dil["cell_indc"][k] = dil["cell_indc"][k - 1]
    print("Finding the cell indices between simulation and reporting grids")
    with alive_bar(len(dil["refxgrid"])) as bar_animation:
        for n, (xcen, zcen) in enumerate(zip(dil["refxgrid"], dil["refzgrid"])):
            bar_animation()
            dil["cell_cent"][n] = np.argmin(
                np.abs(dil["simxcent"] - xcen) + np.abs(dil["simzcent"] - zcen)
            )


def handle_yaxis_mapping_extensive(dig, dil):
    """
    Extend the indices accounting for the y direction (extensive quantities)

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary

    Returns:
        dil (dict): Modified local dictionary

    """
    simyvert = [0]
    for ycent in dil["simycent"]:
        simyvert.append(simyvert[-1] + 2 * (ycent - simyvert[-1]))
    weights, indy, ind = [], [], 0
    for i, (y_i, y_f) in enumerate(zip(simyvert[:-1], simyvert[1:])):
        if dil["refyvert"][ind + 1] <= y_i:
            ind += 1
        if dil["refyvert"][ind] <= y_i and y_f <= dil["refyvert"][ind + 1]:
            indy.append(ind)
            weights.append([1.0])
        else:
            indy.append(ind)
            weights.append([])
            weights[-1].append((dil["refyvert"][ind + 1] - y_i) / (y_f - y_i))
            weights[-1].append((y_f - dil["refyvert"][ind + 1]) / (y_f - y_i))
            ind += 1
    wei_inds = [[] for _ in range(dig["nocellst"])]
    for indz in range(dig["gxyz"][2]):
        i_i = dig["gxyz"][0] * (dig["gxyz"][2] - indz - 1)
        iii = dig["gxyz"][0] * dig["gxyz"][1] * (dig["gxyz"][2] - 1 - indz)
        maps = [
            [
                [
                    col[0]
                    + int(np.floor(col[0] / dig["nxyz"][0]))
                    * dig["nxyz"][0]
                    * (dig["nxyz"][1] - 1),
                    col[1] * weights[0][0],
                ]
                for col in row
            ]
            for row in dil["cell_ind"][i_i : i_i + dig["gxyz"][0]]
        ]
        wei_inds[iii : iii + dig["gxyz"][0]] = maps
        for j, iy in enumerate(indy[1:]):
            wei_inds[
                iii + dig["gxyz"][0] * (j + 1) : iii + dig["gxyz"][0] * (j + 2)
            ] = [
                [
                    [col[0] + (iy * dig["nxyz"][0]), col[1] * weights[j + 1][0]]
                    for col in row
                ]
                for row in maps
            ]
            for weight in weights[j + 1][1:]:
                for i, val in enumerate(
                    wei_inds[
                        iii + dig["gxyz"][0] * (j + 1) : iii + dig["gxyz"][0] * (j + 2)
                    ]
                ):
                    if wei_inds[
                        iii + dig["gxyz"][0] * (j + 1) : iii + dig["gxyz"][0] * (j + 2)
                    ][i]:
                        wei_inds[
                            iii
                            + dig["gxyz"][0] * (j + 1) : iii
                            + dig["gxyz"][0] * (j + 2)
                        ][i].append(
                            [
                                (
                                    val[0][0] + dig["nxyz"][0]
                                    if val[0][0] + dig["nxyz"][0] < dig["nocellsr"]
                                    else val[0][0]
                                ),
                                weight,
                            ]
                        )
    dil["cell_ind"] = wei_inds


def handle_yaxis_mapping_intensive(dig, dil):
    """
    Extend the indices accounting for the y direction

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary

    Returns:
        dil (dict): Modified local dictionary

    """
    indy = np.array(
        [np.argmin(np.abs(dil["simycent"] - y_c)) for y_c in dil["refycent"]]
    )
    tmp_inds = np.zeros(dig["nocellsr"], dtype=int)
    mults = np.zeros(dig["nxyz"][0], dtype=int)
    dil["cell_cent"] = 1.0 * np.array(dil["cell_cent"])
    for indz in range(dig["nxyz"][2]):
        i_i = dig["nxyz"][0] * (dig["nxyz"][2] - indz - 1)
        iii = dig["nxyz"][0] * dig["nxyz"][1] * (dig["nxyz"][2] - 1 - indz)
        if indz != 0:
            mults = np.floor(
                dil["cell_cent"][i_i : i_i + dig["nxyz"][0]] / dig["gxyz"][0]
            )
        values = dil["cell_cent"][i_i : i_i + dig["nxyz"][0]] + mults * (
            dig["gxyz"][0]
        ) * (dig["gxyz"][1] - 1)
        tmp_inds[iii : iii + dig["nxyz"][0]] = values
        for j, iy in enumerate(indy[1:]):
            tmp_inds[
                iii + dig["nxyz"][0] * (j + 1) : iii + dig["nxyz"][0] * (j + 2)
            ] = (iy * dig["gxyz"][0]) + values
    dil["cell_cent"] = tmp_inds


def handle_inactive_mapping(dig, dil):
    """
    Set to inf the inactive grid centers in the reporting grid

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary

    Returns:
        dil (dict): Modified local dictionary

    """
    for i in dig["actind"]:
        for mask in dil["cell_ind"][i]:
            dig["actindr"].append(mask[0])
    dig["actindr"] = list(dict.fromkeys(dig["actindr"]))
    allc = np.linspace(0, dig["nocellsr"] - 1, dig["nocellsr"], dtype=int)
    dig["actindr"] = np.delete(allc, dig["actindr"])


def handle_performance_spatial(dig, dil):
    """
    Create the performance spatial maps

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary

    Returns:
        None

    """
    dil["counter"] = 0.0 * np.ones(dig["nocellsr"])
    dil["pv"] = 0.0 * np.ones(dig["nocellsr"])
    dil["pv"][dig["actindr"]] = 1.0
    static_map_to_report_grid_performance_spatial(dig, dil)
    names = ["co2mn", "h2omn", "co2mb", "h2omb"]
    for i, rst in enumerate(dil["rstno"]):
        if i + 1 < dil["nrstno"]:
            print(
                f"Processing performance spatial {i+1} out of {dil['nrstno']}", end="\r"
            )
        else:
            print(f"Processing performance spatial {i+1} out of {dil['nrstno']}")
        for name in names:
            dil[f"{name}_array"] = np.zeros(dig["nocellst"])
        t_n = rst + dig["no_skip_rst"]
        if t_n > 0:
            # RESIDUAL not included in the SOLUTION deck section (substract 1)
            generate_arrays_performance_spatial(dig, dil, t_n - 1)
        map_to_report_grid_performance_spatial(dig, dil, names, dil["latest_dts"][i])
        write_dense_data_performance_spatial(dig, dil, i)


def static_map_to_report_grid_performance_spatial(dig, dil):
    """
    Map the no dynamic quantities to the reporting grid

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary

    Returns:
        dil (dict): Modified local dictionary

    """
    dil["latest_dts"], infotimes, tsteps = [], [], []
    with open(
        f"{dig['flowf']}/{dig['path'].split('/')[-1].upper()}.INFOSTEP",
        "r",
        encoding="utf8",
    ) as file:
        for j, row in enumerate(csv.reader(file)):
            if j > 0:
                infotimes.append(86400.0 * float((row[0].strip()).split()[0]))
                tsteps.append(86400.0 * float((row[0].strip()).split()[1]))
    infotimes = np.array(infotimes)
    for time in dig["dense_t"][:-1]:
        ind = np.argmin(np.abs(infotimes - (time + dig["time_initial"])))
        if ind > 0:
            dil["latest_dts"].append(tsteps[ind - 1])
        else:
            dil["latest_dts"].append(0.0)
    dil["latest_dts"].append(tsteps[-1])
    for name in ["cvol", "arat"]:
        dil[f"{name}_array"] = np.zeros(dig["nocellst"])
        dil[f"{name}_refg"] = np.zeros(dig["nocellsr"])
    dil["cvol_array"][dig["actind"]] = np.divide(
        dig["porva"], np.array(dig["init"]["PORO"])
    )
    if dig["case"] != "spe11c":
        dil["arat_array"][dig["actind"]] = np.divide(
            np.array(dig["init"]["DZ"]), np.array(dig["init"]["DX"])
        )
    else:
        dil["arat_array"][dig["actind"]] = np.divide(
            np.array(dig["init"]["DZ"]),
            (np.array(dig["init"]["DX"]) ** 2 + np.array(dig["init"]["DY"]) ** 2)
            ** 0.5,
        )
    for i in dig["actind"]:
        for mask in dil["cell_ind"][i]:
            dil["cvol_refg"][mask[0]] += dil["cvol_array"][i]
            dil["arat_refg"][mask[0]] += dil["arat_array"][i]
            dil["counter"][mask[0]] += 1.0
            dil["pv"][mask[0]] += dig["porv"][i]
    inds = dil["counter"] > 0.0
    dil["cvol_refg"][inds] = np.divide(dil["cvol_refg"][inds], dil["counter"][inds])
    dil["arat_refg"][inds] = np.divide(dil["arat_refg"][inds], dil["counter"][inds])
    for name in ["cvol", "arat"]:
        dil[f"{name}_refg"][dil[f"{name}_refg"] < 1e-12] = np.nan
    dil["ei"] = dil["pv"] > 0.0


def generate_arrays_performance_spatial(dig, dil, t_n):
    """
    Arrays for the performance spatial data

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary\n
        t_n (int): Index for the number of restart file

    Returns:
        dil (dict): Modified local dictionary

    """
    dil["co2mb_array"][dig["actind"]] = np.array(dig["unrst"]["RES_GAS", t_n + 1])
    if dig["unrst"].count("RES_WAT", t_n + 1):
        dil["h2omb_array"][dig["actind"]] = np.array(dig["unrst"]["RES_WAT", t_n + 1])
    else:
        dil["h2omb_array"][dig["actind"]] = np.array(dig["unrst"]["RES_OIL", t_n + 1])
    dil["co2mn_array"][dig["actind"]] = np.divide(
        np.abs(dil["co2mb_array"][dig["actind"]]), dig["porva"]
    )
    dil["h2omn_array"][dig["actind"]] = np.divide(
        np.abs(dil["h2omb_array"][dig["actind"]]), dig["porva"]
    )


def map_to_report_grid_performance_spatial(dig, dil, names, d_t):
    """
    Map the simulation grid to the reporting grid

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary\n
        names (list): Strings with the quantities for the spatial maps\n
        d_t (float): Time step size

    Returns:
        dil (dict): Modified local dictionary

    """
    for name in names:
        dil[f"{name}_refg"] = np.zeros(dig["nocellsr"])
        dil[f"{name}_refg"][dig["actindr"]] = np.nan
    for i in dig["actind"]:
        for mask in dil["cell_ind"][i]:
            dil["co2mn_refg"][mask[0]] = max(
                dil["co2mn_refg"][mask[0]], dil["co2mn_array"][i] * mask[1]
            )
            dil["h2omn_refg"][mask[0]] = max(
                dil["h2omn_refg"][mask[0]], dil["h2omn_array"][i] * mask[1]
            )
            dil["co2mb_refg"][mask[0]] += dil["co2mb_array"][i] * mask[1]
            dil["h2omb_refg"][mask[0]] += dil["h2omb_array"][i] * mask[1]
    dil["co2mn_refg"] *= d_t
    dil["h2omn_refg"] *= d_t
    dil["co2mb_refg"][dil["ei"]] = d_t * np.divide(
        dil["co2mb_refg"][dil["ei"]], dil["pv"][dil["ei"]]
    )
    dil["h2omb_refg"][dil["ei"]] = d_t * np.divide(
        dil["h2omb_refg"][dil["ei"]], dil["pv"][dil["ei"]]
    )


def write_dense_data_performance_spatial(dig, dil, i):
    """
    Write the dense performance data

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary\n
        i (int): Number of csv file

    Returns:
        None

    """
    if dig["case"] == "spe11a":
        name_t = f"{round(dig['dense_t'][i]/3600)}h"
    else:
        name_t = f"{round(dig['dense_t'][i]/SECONDS_IN_YEAR)}y"
    if dig["case"] != "spe11c":
        text = [
            "# x [m], z [m], cvol [m^3], arat [-], CO2 max_norm_res [-], "
            + "H2O max_norm_res [-], CO2 mb_error [-], H2O mb_error [-], "
            + "post_est [-]"
        ]
    else:
        text = [
            "# x [m], y [m], z [m], cvol [m^3], arat [-], CO2 max_norm_res [-], "
            + "H2O max_norm_res [-], CO2 mb_error [-], H2O mb_error [-], "
            + "post_est [-]"
        ]
    idz = 0
    for zcord in dil["refzcent"]:
        idxy = 0
        for ycord in dil["refycent"]:
            for xcord in dil["refxcent"]:
                idc = -dig["nxyz"][0] * dig["nxyz"][1] * (dig["nxyz"][2] - idz) + idxy
                if dig["case"] != "spe11c":
                    if np.isnan(dil["cvol_refg"][idc]):
                        text.append(
                            f"{xcord:.3e}, {zcord:.3e}, n/a, n/a, n/a, n/a, n/a, "
                            + "n/a, n/a"
                        )
                    else:
                        text.append(
                            f"{xcord:.3e}, {zcord:.3e}, "
                            + f"{dil['cvol_refg'][idc] :.3e}, "
                            + f"{dil['arat_refg'][idc] :.3e}, "
                            + f"{dil['co2mn_refg'][idc] :.3e}, "
                            + f"{dil['h2omn_refg'][idc] :.3e}, "
                            + f"{dil['co2mb_refg'][idc] :.3e}, "
                            + f"{dil['h2omb_refg'][idc] :.3e}, "
                            + "n/a"
                        )
                else:
                    if np.isnan(dil["cvol_refg"][idc]):
                        text.append(
                            f"{xcord:.3e}, {ycord:.3e}, {zcord:.3e}, n/a, n/a, n/a, "
                            + "n/a, n/a, n/a, n/a"
                        )
                    else:
                        text.append(
                            f"{xcord:.3e}, {ycord:.3e}, {zcord:.3e}, "
                            + f"{dil['cvol_refg'][idc] :.3e}, "
                            + f"{dil['arat_refg'][idc] :.3e}, "
                            + f"{dil['co2mn_refg'][idc] :.3e}, "
                            + f"{dil['h2omn_refg'][idc] :.3e}, "
                            + f"{dil['co2mb_refg'][idc] :.3e}, "
                            + f"{dil['h2omb_refg'][idc] :.3e}, "
                            + "n/a"
                        )
                idxy += 1
        idz += 1
    with open(
        f"{dig['where']}/{dig['case']}_performance_spatial_map_{name_t}.csv",
        "w",
        encoding="utf8",
    ) as file:
        file.write("\n".join(text))


def generate_arrays(dig, dil, names, t_n):
    """
    Arrays for the dense data

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary\n
        names (list): Strings with the quantities for the spatial maps\n
        t_n (int): Index for the number of restart file

    Returns:
        dil (dict): Modified local dictionary

    """
    for name in names[:-1]:
        dil[f"{name}_array"] = np.zeros(dig["nocellst"])
        dil[f"{name}_refg"] = np.zeros(dig["nocellsr"])
        if dig["case"] == "spe11a" or (dig["lower"] and not dig["cp"]):
            dil[f"{name}_array"] = np.empty(dig["nocellst"]) * np.nan
    dil["tco2_array"] = np.zeros(dig["nocellst"])
    dil["tco2_refg"] = np.zeros(dig["nocellsr"])
    dil["tco2_refg"][dig["actindr"]] = np.nan
    sgas = np.abs(np.array(dig["unrst"]["SGAS", t_n]))
    rhog = np.array(dig["unrst"]["GAS_DEN", t_n])
    pres = np.array(dig["unrst"]["PRESSURE", t_n]) - np.array(dig["unrst"]["PCGW", t_n])
    rhow = np.array(dig["unrst"]["WAT_DEN", t_n])
    rvv, rss = 0.0 * sgas, 0.0 * sgas
    if not dig["immiscible"]:
        rss = np.array(dig["unrst"]["RSW", t_n])
        if dig["unrst"].count("RVW", t_n):
            rvv = np.array(dig["unrst"]["RVW", t_n])
        if dig["case"] != "spe11a":
            dil["temp_array"][dig["actind"]] = np.array(dig["unrst"]["TEMP", t_n])
        x_l_co2 = np.divide(rss, rss + WAT_DEN_REF / GAS_DEN_REF)
        x_g_h2o = np.divide(rvv, rvv + GAS_DEN_REF / WAT_DEN_REF)
        co2_g = (1 - x_g_h2o) * sgas * rhog * dig["porva"]
        co2_d = x_l_co2 * (1 - sgas) * rhow * dig["porva"]
    else:
        x_l_co2, x_g_h2o, co2_d = 0, 0, 0
        co2_g = sgas * rhog * dig["porva"]
    dil["pressure_array"][dig["actind"]] = 1e5 * pres
    dil["sgas_array"][dig["actind"]] = sgas * (sgas > SGAS_THR)
    dil["gden_array"][dig["actind"]] = rhog * (sgas > SGAS_THR)
    dil["wden_array"][dig["actind"]] = rhow
    dil["xco2_array"][dig["actind"]] = x_l_co2
    dil["xh20_array"][dig["actind"]] = x_g_h2o * (sgas > SGAS_THR)
    dil["tco2_array"][dig["actind"]] = co2_d + co2_g
    if dig["lower"] and dig["cp"]:
        for name in ["pressure", "sgas", "gden", "wden", "xco2", "xh20", "temp"]:
            if f"{name}_array" in dil:
                dil[f"{name}_array"] = np.insert(
                    dil[f"{name}_array"],
                    0,
                    np.full(shape=dig["gxyz"][0] * dig["gxyz"][1], fill_value=np.nan),
                )


def map_to_report_grid(dig, dil, names):
    """
    Map the simulation grid to the reporting grid

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary\n
        names (list): Strings with the quantities for the spatial maps

    Returns:
        dil (dict): Modified local dictionary

    """
    for i in dig["actind"]:
        for mask in dil["cell_ind"][i]:
            dil["tco2_refg"][mask[0]] += dil["tco2_array"][i] * mask[1]
    for i, ind in enumerate(dil["cell_cent"]):
        for name in names[:-1]:
            dil[f"{name}_refg"][i] = dil[f"{name}_array"][int(ind)]


def write_dense_data(dig, dil, i):
    """
    Map the quantities to the cells

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary\n
        i (int): Number of csv file

    Returns:
        None

    """
    name_t, text = get_header(dig, i)
    idz = 0
    if dig["lower"]:
        dil["tco2_refg"][np.isnan(dil["sgas_refg"])] = np.nan
    for zcord in dil["refzcent"]:
        idxy = 0
        for ycord in dil["refycent"]:
            for xcord in dil["refxcent"]:
                idc = -dig["nxyz"][0] * dig["nxyz"][1] * (dig["nxyz"][2] - idz) + idxy
                co2, xco2, xh20, temp = "n/a", "n/a", "n/a", "n/a"
                if not np.isnan(dil["tco2_refg"][idc]):
                    co2 = f"{dil['tco2_refg'][idc] :.3e}"
                if not dig["immiscible"]:
                    xco2 = f"{dil['xco2_refg'][idc] :.3e}"
                    xh20 = f"{dil['xh20_refg'][idc] :.3e}"
                    if dig["case"] != "spe11a":
                        temp = f"{dil['temp_refg'][idc] :.3e}"
                if dig["case"] == "spe11a":
                    if np.isnan(dil["pressure_refg"][idc]):
                        text.append(
                            f"{xcord:.3e}, {zcord:.3e}, n/a, n/a, n/a, n/a, n/a, "
                            + f"n/a, {co2}"
                        )
                    else:
                        text.append(
                            f"{xcord:.3e}, {zcord:.3e}, "
                            + f"{dil['pressure_refg'][idc] :.3e}, "
                            + f"{dil['sgas_refg'][idc] :.3e}, "
                            + f"{xco2}, "
                            + f"{xh20}, "
                            + f"{dil['gden_refg'][idc] :.3e}, "
                            + f"{dil['wden_refg'][idc] :.3e}, "
                            + f"{co2}"
                        )
                elif dig["case"] == "spe11b":
                    if np.isnan(dil["pressure_refg"][idc]):
                        text.append(
                            f"{xcord:.3e}, {zcord:.3e}, "
                            + f"n/a, n/a, n/a, n/a, n/a, n/a, {co2}, n/a"
                        )
                    else:
                        text.append(
                            f"{xcord:.3e}, {zcord:.3e}, "
                            + f"{dil['pressure_refg'][idc] :.3e}, "
                            + f"{dil['sgas_refg'][idc] :.3e}, "
                            + f"{xco2}, "
                            + f"{xh20}, "
                            + f"{dil['gden_refg'][idc] :.3e}, "
                            + f"{dil['wden_refg'][idc] :.3e}, "
                            + f"{co2}, "
                            + f"{temp}"
                        )
                else:
                    if np.isnan(dil["pressure_refg"][idc]):
                        text.append(
                            f"{xcord:.3e}, {ycord:.3e}, {zcord:.3e}, "
                            + f"n/a, n/a, n/a, n/a, n/a, n/a, {co2}, n/a"
                        )
                    else:
                        text.append(
                            f"{xcord:.3e}, {ycord:.3e}, {zcord:.3e}, "
                            + f"{dil['pressure_refg'][idc] :.3e}, "
                            + f"{dil['sgas_refg'][idc] :.3e}, "
                            + f"{xco2}, "
                            + f"{xh20}, "
                            + f"{dil['gden_refg'][idc] :.3e}, "
                            + f"{dil['wden_refg'][idc] :.3e}, "
                            + f"{co2}, "
                            + f"{temp}"
                        )
                idxy += 1
        idz += 1
    with open(
        f"{dig['where']}/{dig['case']}_spatial_map_{name_t}.csv", "w", encoding="utf8"
    ) as file:
        file.write("\n".join(text))


def get_header(dig, i):
    """
    Get the right csv file header

    Args:
        dig (dict): Global dictionary\n
        i (int): Number of csv file

    Returns:
        name_t (str): Name of the csv file\n
        text (str): Header for the csv file

    """
    if dig["case"] == "spe11a":
        name_t = f"{round(dig['dense_t'][i]/3600)}h"
        text = [
            "# x [m], z [m], pressure [Pa], gas saturation [-], "
            + "mass fraction of CO2 in liquid [-], mass fraction of H20 in vapor [-], "
            + "phase mass density gas [kg/m3], phase mass density water [kg/m3], "
            + "total mass CO2 [kg]"
        ]
    elif dig["case"] == "spe11b":
        name_t = f"{round(dig['dense_t'][i]/SECONDS_IN_YEAR)}y"
        text = [
            "# x [m], z [m], pressure [Pa], gas saturation [-], "
            + "mass fraction of CO2 in liquid [-], mass fraction of H20 in vapor [-], "
            + "phase mass density gas [kg/m3], phase mass density water [kg/m3], "
            + "total mass CO2 [kg], temperature [C]"
        ]
    else:
        name_t = f"{round(dig['dense_t'][i]/SECONDS_IN_YEAR)}y"
        text = [
            "# x [m], y [m], z [m], pressure [Pa], gas saturation [-], "
            + "mass fraction of CO2 in liquid [-], mass fraction of H20 in vapor [-], "
            + "phase mass density gas [kg/m3], phase mass density water [kg/m3], "
            + "total mass CO2 [kg], temperature [C]"
        ]
    return name_t, text


if __name__ == "__main__":
    main()
