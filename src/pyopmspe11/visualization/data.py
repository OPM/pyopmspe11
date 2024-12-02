# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT
# pylint: disable=C0302, R0912, R0914, R0801

""""
Script to write the benchmark data
"""

import os
import argparse
import warnings
import csv
from io import StringIO
from shapely.geometry import Polygon
from rtree import index
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from resdata.grid import Grid
from resdata.resfile import ResdataFile
from resdata.summary import Summary

try:
    from opm.io.ecl import EclFile as OpmFile
    from opm.io.ecl import EGrid as OpmGrid
    from opm.io.ecl import ERst as OpmRestart
    from opm.io.ecl import ESmry as OpmSummary
except ImportError:
    pass

GAS_DEN_REF = 1.86843
WAT_DEN_REF = 998.108
SECONDS_IN_YEAR = 31536000
KMOL_TO_KG = 1e3 * 0.044
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
        "-u",
        "--use",
        default="resdata",
        help="Using the 'resdata' or python package (resdata by default).",
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
    dig["exe"] = os.getcwd()
    dig["where"] = f"{dig['exe']}/{dig['path']}/data"
    dig["use"] = cmdargs["use"].strip()
    dig["nxyz"] = np.genfromtxt(
        StringIO(cmdargs["resolution"]), delimiter=",", dtype=int
    )
    dig["sim"] = "./" + dig["path"] + "/flow/" + f"{dig['path'].upper()}"
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
    read_times(dig)
    if dig["use"] == "opm":
        read_opm(dig)
    else:
        read_resdata(dig)
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


def read_times(dig):
    """
    Get the time for injection and restart number

    Args:
        dig (dict): Global dictionary

    Returns:
        dig (dict): Modified global dictionary

    """
    with open(f"{dig['path']}/deck/dt.txt", "r", encoding="utf8") as file:
        for i, value in enumerate(csv.reader(file)):
            if i == 0:
                dig["time_initial"] = float(value[0])
            elif i == 1:
                dig["no_skip_rst"] = int(value[0])
            else:
                dig["times"] = list(
                    np.genfromtxt(StringIO(value[0]), delimiter=" ", dtype=float)
                )


def read_resdata(dig):
    """
    Read the simulation files using resdata

    Args:
        dig (dict): Global dictionary

    Returns:
        dig (dict): Modified global dictionary

    """
    dig["unrst"] = ResdataFile(f"{dig['sim']}.UNRST")
    dig["init"] = ResdataFile(f"{dig['sim']}.INIT")
    dig["egrid"] = Grid(f"{dig['sim']}.EGRID")
    dig["smspec"] = Summary(f"{dig['sim']}.SMSPEC")
    if dig["unrst"].has_kw("WAT_DEN"):
        dig["watDen"], dig["r_s"], dig["r_v"] = "wat_den", "rsw", "rvw"
        dig["bpr"] = "BWPR"
    else:
        dig["watDen"], dig["r_s"], dig["r_v"] = "oil_den", "rs", "rv"
        dig["bpr"] = "BPR"
    dig["porv"] = np.array(dig["init"].iget_kw("PORV")[0])
    dig["actnum"] = list(dig["egrid"].export_actnum())
    dig["actind"] = list(i for i, act in enumerate(dig["actnum"]) if act == 1)
    dig["porva"] = np.array(
        [porv for (porv, act) in zip(dig["porv"], dig["actnum"]) if act == 1]
    )
    dig["nocellst"], dig["nocellsa"] = len(dig["actnum"]), sum(dig["actnum"])
    dig["norst"] = dig["unrst"].num_report_steps()
    dig["times_summary"] = [0.0]
    dig["times_summary"] += list(
        86400.0 * dig["smspec"]["TIME"].values - dig["time_initial"]
    )
    dig["gxyz"] = [
        dig["egrid"].nx,
        dig["egrid"].ny,
        dig["egrid"].nz,
    ]
    dig["noxz"] = dig["egrid"].nx * dig["egrid"].nz


def read_opm(dig):
    """
    Read the simulation files using OPM

    Args:
        dig (dict): Global dictionary

    Returns:
        dig (dict): Modified global dictionary

    """
    dig["unrst"] = OpmRestart(f"{dig['sim']}.UNRST")
    dig["init"] = OpmFile(f"{dig['sim']}.INIT")
    dig["egrid"] = OpmGrid(f"{dig['sim']}.EGRID")
    dig["smspec"] = OpmSummary(f"{dig['sim']}.SMSPEC")
    opm_files(dig)
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


def opm_files(dig):
    """
    Read some of the data from the simulation output files

    Args:
        dig (dict): Global dictionary

    Returns:
        dig (dict): Modified global dictionary

    """
    dig["norst"] = len(dig["unrst"].report_steps)
    if dig["unrst"].count("WAT_DEN", 0):
        dig["watDen"], dig["r_s"], dig["r_v"] = "wat_den", "rsw", "rvw"
        dig["bpr"] = "BWPR"
    else:
        dig["watDen"], dig["r_s"], dig["r_v"] = "oil_den", "rs", "rv"
        dig["bpr"] = "BPR"
    dig["porv"] = np.array(dig["init"]["PORV"])


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
        f"{dig['path']}/flow/{dig['path'].upper()}.INFOSTEP", "r", encoding="utf8"
    ) as file:
        for j, row in enumerate(csv.reader(file)):
            if j > 0:
                if float((row[0].strip()).split()[0]) >= (
                    (dig["time_initial"] - dig["sparse_t"]) / 86400.0
                ):
                    dil["infosteps"].append(
                        [float(column) for column in (row[0].strip()).split()]
                    )
    infotimes = [
        infostep[0] * 86400 - dig["time_initial"] for infostep in dil["infosteps"]
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
    for j in range(max(dil["detail_info"]) + 1):
        ind = j == dil["detail_info"]
        dil["times_det"].append(max(infotimes[ind]))
    dil["times_det"] = np.array(dil["times_det"])
    dil["fsteps"] = np.array(
        [1.0 * (infostep[11] == 0) for infostep in dil["infosteps"]]
    )
    dil["nress"] = np.array([infostep[8] for infostep in dil["infosteps"]])
    dil["tlinsols"] = np.array([infostep[4] for infostep in dil["infosteps"]])
    dil["liniters"] = np.array([infostep[10] for infostep in dil["infosteps"]])
    dil["nliters"] = np.array([infostep[9] for infostep in dil["infosteps"]])
    dil["tsteps"] = np.array(
        [86400 * infostep[1] * infostep[11] for infostep in dil["infosteps"]]
    )
    dil["alltsteps"] = np.array([86400 * infostep[1] for infostep in dil["infosteps"]])
    if dig["use"] == "opm":
        tcpu = dig["smspec"]["TCPU"]
        fgip = GAS_DEN_REF * dig["smspec"]["FGIP"]
        times = 86400.0 * dig["smspec"]["TIME"] - dig["time_initial"]
    else:
        tcpu = dig["smspec"]["TCPU"].values
        fgip = GAS_DEN_REF * dig["smspec"]["FGIP"].values
        times = 86400.0 * dig["smspec"]["TIME"].values - dig["time_initial"]
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
        fgip = np.insert(fgip, 0, 0)
    interp_fgip = interp1d(
        times,
        fgip,
        fill_value="extrapolate",
    )
    write_performance(dig, dil, interp_fgip, tcpu, infotimes)


def write_performance(dig, dil, interp_fgip, tcpu, infotimes):
    """
    Write the performance data

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary\n
        interp_fgip (object): Interpolator (time) for the CO2 mass\n
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
        if sum(tcpu[itd]) == 0:
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
            + f"{interp_fgip(time):.3e}, "
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
    for j in range(max(dil["detail_info"]) + 1):
        ind = j == dil["detail_info"]
        time = max(infotimes[ind])
        if time >= 0:
            dil["text"].append(
                f"{time:.3e}, "
                + f"{max(dil['tsteps'][ind]):.3e}, "
                + f"{np.sum(dil['fsteps'][ind]):.3e}, "
                + f"{interp_fgip(time):.3e}, "
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
        if key[0 : len(dig["bpr"])] == dig["bpr"] and "," in key[len(dig["bpr"]) + 1 :]:
            names.append(key)
            i_jk.append(
                np.genfromtxt(
                    StringIO(key[len(dig["bpr"]) + 1 :]), delimiter=",", dtype=int
                )[0]
            )
            ind += 1
            if ind == 2:
                break
    sort = sorted(range(len(i_jk)), key=i_jk.__getitem__)
    if dig["use"] == "opm":
        pop1 = dig["unrst"]["PRESSURE", 0][dil["fipnum"].index(8)]
        pop2 = dig["unrst"]["PRESSURE", 0][dil["fipnum"].index(9)]
        if dig["unrst"].count("PCGW", 0):
            pop1 -= dig["unrst"]["PCGW", 0][dil["fipnum"].index(8)]
            pop2 -= dig["unrst"]["PCGW", 0][dil["fipnum"].index(9)]
        dil["pop1"] = [pop1 * 1.0e5] + list(dig["smspec"][names[sort[0]]] * 1.0e5)  # Pa
        dil["pop2"] = [pop2 * 1.0e5] + list(dig["smspec"][names[sort[1]]] * 1.0e5)  # Pa
        for i in dil["fip_diss_a"]:
            dil["moba"] += dig["smspec"][f"RGKDM:{i}"] * KMOL_TO_KG
            dil["imma"] += dig["smspec"][f"RGKDI:{i}"] * KMOL_TO_KG
            dil["dissa"] += dig["smspec"][f"RWCD:{i}"] * KMOL_TO_KG
        for i in dil["fip_seal_a"]:
            dil["seala"] += (
                dig["smspec"][f"RWCD:{i}"]
                + dig["smspec"][f"RGKDM:{i}"]
                + dig["smspec"][f"RGKDI:{i}"]
            ) * KMOL_TO_KG
        for i in dil["fip_diss_b"]:
            dil["mobb"] += dig["smspec"][f"RGKDM:{i}"] * KMOL_TO_KG
            dil["immb"] += dig["smspec"][f"RGKDI:{i}"] * KMOL_TO_KG
            dil["dissb"] += dig["smspec"][f"RWCD:{i}"] * KMOL_TO_KG
        for i in dil["fip_seal_b"]:
            dil["sealb"] += (
                dig["smspec"][f"RWCD:{i}"]
                + dig["smspec"][f"RGKDM:{i}"]
                + dig["smspec"][f"RGKDI:{i}"]
            ) * KMOL_TO_KG
        dil["sealt"] = dil["seala"] + dil["sealb"]
        for name in ["RWCD", "RGKDM", "RGKDI"]:
            dil["sealt"] += (
                dig["smspec"][f"{name}:7"] + dig["smspec"][f"{name}:9"]
            ) * KMOL_TO_KG
        if dig["case"] != "spe11a":
            sealbound = (
                dig["smspec"]["RWCD:10"]
                + dig["smspec"]["RGKDM:10"]
                + dig["smspec"]["RGKDI:10"]
            ) * KMOL_TO_KG
            dil["sealt"] += sealbound
            dil["boundtot"] = sealbound
            for i in dil["fip_bound_t"]:
                dil["boundtot"] += (
                    dig["smspec"][f"RWCD:{i}"]
                    + dig["smspec"][f"RGKDM:{i}"]
                    + dig["smspec"][f"RGKDI:{i}"]
                ) * KMOL_TO_KG
    else:
        resdata_summary(dig, dil, names, sort)


def resdata_summary(dig, dil, names, sort):
    """
    Read the summary arrays using resdata

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary\n
        names (list): Strings with the sensors ijk locations\n
        sort (list): Integers with the right order for the sensors

    Returns:
        dil (dict): Modified local dictionary

    """
    pop1 = dig["unrst"]["PRESSURE"][0][dil["fipnum"].index(8)]
    pop2 = dig["unrst"]["PRESSURE"][0][dil["fipnum"].index(9)]
    if dig["unrst"].has_kw("PCGW"):
        pop1 -= dig["unrst"]["PCGW"][0][dil["fipnum"].index(8)]
        pop2 -= dig["unrst"]["PCGW"][0][dil["fipnum"].index(9)]
    dil["pop1"] = [pop1 * 1.0e5] + list(
        dig["smspec"][names[sort[0]]].values * 1.0e5
    )  # Pa
    dil["pop2"] = [pop2 * 1.0e5] + list(
        dig["smspec"][names[sort[1]]].values * 1.0e5
    )  # Pa
    for i in dil["fip_diss_a"]:
        dil["moba"] += dig["smspec"][f"RGKDM:{i}"].values * KMOL_TO_KG
        dil["imma"] += dig["smspec"][f"RGKDI:{i}"].values * KMOL_TO_KG
        dil["dissa"] += dig["smspec"][f"RWCD:{i}"].values * KMOL_TO_KG
    for i in dil["fip_seal_a"]:
        dil["seala"] += (
            dig["smspec"][f"RWCD:{i}"].values
            + dig["smspec"][f"RGKDM:{i}"].values
            + dig["smspec"][f"RGKDI:{i}"].values
        ) * KMOL_TO_KG
    for i in dil["fip_diss_b"]:
        dil["mobb"] += dig["smspec"][f"RGKDM:{i}"].values * KMOL_TO_KG
        dil["immb"] += dig["smspec"][f"RGKDI:{i}"].values * KMOL_TO_KG
        dil["dissb"] += dig["smspec"][f"RWCD:{i}"].values * KMOL_TO_KG
    for i in dil["fip_seal_b"]:
        dil["sealb"] += (
            dig["smspec"][f"RWCD:{i}"].values
            + dig["smspec"][f"RGKDM:{i}"].values
            + dig["smspec"][f"RGKDI:{i}"].values
        ) * KMOL_TO_KG
    dil["sealt"] = dil["seala"] + dil["sealb"]
    for name in ["RWCD", "RGKDM", "RGKDI"]:
        dil["sealt"] += (
            dig["smspec"][f"{name}:7"].values + dig["smspec"][f"{name}:9"].values
        ) * KMOL_TO_KG
    if dig["case"] != "spe11a":
        sealbound = (
            dig["smspec"]["RWCD:10"].values
            + dig["smspec"]["RGKDM:10"].values
            + dig["smspec"]["RGKDI:10"].values
        ) * KMOL_TO_KG
        dil["sealt"] += sealbound
        dil["boundtot"] = sealbound
        for i in dil["fip_bound_t"]:
            dil["boundtot"] += (
                dig["smspec"][f"RWCD:{i}"].values
                + dig["smspec"][f"RGKDM:{i}"].values
                + dig["smspec"][f"RGKDI:{i}"].values
            ) * KMOL_TO_KG


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
    if dig["use"] == "opm":
        dil["fipnum"] = list(dig["init"]["FIPNUM"])
        for name in ["dx", "dy", "dz"]:
            dil[f"{name}"] = np.array(dig["init"][name.upper()])
    else:
        dil["fipnum"] = list(dig["init"].iget_kw("FIPNUM")[0])
        for name in ["dx", "dy", "dz"]:
            dil[f"{name}"] = np.array(dig["init"].iget_kw(name.upper())[0])
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
    compute_m_c(dig, dil)
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
    dil["fip_diss_a"] = [2, 4, 5, 8]
    dil["fip_seal_a"] = [5, 8]
    dil["fip_diss_b"] = [3, 6]
    dil["fip_seal_b"] = [6]
    if dig["case"] != "spe11a":
        dil["fip_bound_t"] = [11]
    if dig["case"] == "spe11c":
        dil["fip_diss_a"] += [13, 14, 17]
        dil["fip_seal_a"] += [14]
        dil["fip_diss_b"] += [15, 16]
        dil["fip_seal_b"] += [16]
        dil["fip_bound_t"] += [13, 14, 15, 16, 17]
        if max(dil["fipnum"]) == 18:
            dil["fip_diss_a"] += [12, 18]
            dil["fip_seal_a"] += [12, 18]
            dil["fip_bound_t"] += [18]


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
    max_xcw(dig, dil)
    for t_n in range(dig["no_skip_rst"] + 1, dig["norst"]):
        if dig["use"] == "opm":
            rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}", t_n])
        else:
            rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}"][t_n])
        dil["xcw"] = np.divide(rss, rss + WAT_DEN_REF / GAS_DEN_REF)
        if dil["xcw_max"] > 0:
            dil["xcw"] /= dil["xcw_max"]
        if dil["xcw_max"] == -1:
            if dig["use"] == "opm":
                rssat = np.array(dig["unrst"][f"{dig['r_s'].upper()}SAT", t_n])
            else:
                rssat = np.array(dig["unrst"][f"{dig['r_s'].upper()}SAT"][t_n])
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
            interp = interp1d(
                dig["times_summary"],
                [0.0] + list(dil[f"{name}"]),
                fill_value="extrapolate",
            )
        dil[f"{name}"] = interp(dil["times_data"])
    text = [
        "# t [s], p1 [Pa], p2 [Pa], mobA [kg], immA [kg], dissA [kg], sealA [kg], "
        + "<same for B>, MC [m^2], sealTot [kg]"
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


def max_xcw(dig, dil):
    """
    Get the maximum CO2 mass fraction in the liquid phase

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary

    Returns:
        dil (dict): Modified local dictionary

    """
    dil["xcw_max"] = 0
    if dig["use"] == "opm":
        if dig["unrst"].count(f"{dig['r_s'].upper()}SAT", 0):
            dil["xcw_max"] = -1
            return
    else:
        if dig["unrst"].has_kw(f"{dig['r_s'].upper()}SAT"):
            dil["xcw_max"] = -1
            return
    for t_n in range(dig["no_skip_rst"], dig["norst"]):
        if dig["use"] == "opm":
            rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}", t_n])
        else:
            rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}"][t_n])
        xcw = np.divide(rss, rss + WAT_DEN_REF / GAS_DEN_REF)
        xcw_max = np.max(xcw[dil["boxc"]])
        dil["xcw_max"] = max(xcw_max, dil["xcw_max"])


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
    with open(f"{dig['path']}/deck/centers.txt", "r", encoding="utf8") as file:
        for j, row in enumerate(csv.reader(file)):
            dil["simxcent"][j] = float(row[0])
            dil["simzcent"][j] = dig["dims"][2] - float(row[2])
            if dil["simzcent"][j] == 0:
                dil["simxcent"][j] = -1e10
                dil["simzcent"][j] = -1e10
    dil["simpoly"] = []
    with open(f"{dig['path']}/deck/corners.txt", "r", encoding="utf8") as file:
        for row in csv.reader(file):
            dil["simpoly"].append(
                Polygon(
                    [
                        (float(row[0]), float(row[1])),
                        (float(row[2]), float(row[3])),
                        (float(row[4]), float(row[5])),
                        (float(row[6]), float(row[7])),
                    ]
                )
            )
    if dig["use"] == "opm":
        dil["satnum"] = list(dig["init"]["SATNUM"])
    else:
        dil["satnum"] = list(dig["init"].iget_kw("SATNUM")[0])


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
    ind, dil["cell_ind"] = 0, [[] for _ in range(dig["noxz"])]
    dil["cell_indc"] = np.zeros(dig["noxz"])
    dil["cell_cent"] = [0 for _ in range(dig["noxzr"])]
    idx = index.Index()
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
    for k, simp in enumerate(dil["simpoly"]):
        ovrl = list(idx.intersection(simp.bounds))
        if simp.area > 0:
            for ind in ovrl:
                area = simp.intersection(dil["refpoly"][ind]).area / simp.area
                if area > 0:
                    dil["cell_ind"][k].append([ind, area])
                    dil["cell_indc"][k] = ind
        else:
            dil["cell_indc"][k] = dil["cell_indc"][k - 1]
    for k, (xcen, zcen) in enumerate(zip(dil["refxgrid"], dil["refzgrid"])):
        dil["cell_cent"][k] = pd.Series(
            np.abs(dil["simxcent"] - xcen) + np.abs(dil["simzcent"] - zcen)
        ).argmin()
    dig["actindr"] = []
    if max(dil["satnum"]) < 7 and dig["case"] == "spe11a":
        handle_inactive_mapping(dig, dil)
    if dig["case"] == "spe11c":
        handle_yaxis_mapping_intensive(dig, dil)
        handle_yaxis_mapping_extensive(dig, dil)
    if dig["mode"] == "all" or dig["mode"][:5] == "dense":
        names = ["pressure", "sgas", "xco2", "xh20", "gden", "wden", "tco2"]
        if dig["case"] != "spe11a":
            names = ["temp"] + names
        for i, rst in enumerate(dil["rstno"]):
            print(f"Processing dense data {i+1} out of {dil['nrstno']}")
            t_n = rst + dig["no_skip_rst"]
            generate_arrays(dig, dil, names, t_n)
            map_to_report_grid(dig, dil, names)
            write_dense_data(dig, dil, i)
    if dig["mode"] in ["all", "performance-spatial", "dense_performance-spatial"]:
        handle_performance_spatial(dig, dil)


def handle_yaxis_mapping_extensive(dig, dil):
    """
    Extend the indices accounting for the y direction (extensive quantities)

    Args:
        dig (dict): Global dictionary\n
        dil (dict): Local dictionary

    Returns:
        dil (dict): Modified local dictionary

    """
    simycent = [0.0] * dig["gxyz"][1]
    with open(f"{dig['path']}/deck/ycenters.txt", "r", encoding="utf8") as file:
        for j, row in enumerate(csv.reader(file)):
            simycent[j] = float(f"{float(row[0]):.1f}")
    simyvert = [0]
    for ycent in simycent:
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
            for i, row in enumerate(dil["cell_ind"][i_i : i_i + dig["gxyz"][0]])
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
    simycent = [0.0] * dig["gxyz"][1]
    with open(f"{dig['path']}/deck/ycenters.txt", "r", encoding="utf8") as file:
        for j, row in enumerate(csv.reader(file)):
            simycent[j] = float(row[0])
    indy = np.array(
        [pd.Series(np.abs(simycent - y_c)).argmin() for y_c in dil["refycent"]]
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
        f"{dig['path']}/flow/{dig['path'].upper()}.INFOSTEP", "r", encoding="utf8"
    ) as file:
        for j, row in enumerate(csv.reader(file)):
            if j > 0:
                infotimes.append(86400.0 * float((row[0].strip()).split()[0]))
                tsteps.append(86400.0 * float((row[0].strip()).split()[1]))
    infotimes = np.array(infotimes)
    for time in dig["dense_t"][:-1]:
        ind = pd.Series(np.abs(infotimes - (time + dig["time_initial"]))).argmin()
        if ind > 0:
            dil["latest_dts"].append(tsteps[ind - 1])
        else:
            dil["latest_dts"].append(0.0)
    dil["latest_dts"].append(tsteps[-1])
    for name in ["cvol", "arat"]:
        dil[f"{name}_array"] = np.zeros(dig["nocellst"])
        dil[f"{name}_refg"] = np.zeros(dig["nocellsr"])
    if dig["use"] == "opm":
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
    else:
        dil["cvol_array"][dig["actind"]] = np.divide(
            dig["porva"], np.array(dig["init"].iget_kw("PORO")[0])
        )
        if dig["case"] != "spe11c":
            dil["arat_array"][dig["actind"]] = np.divide(
                np.array(dig["init"].iget_kw("DZ")[0]),
                np.array(dig["init"].iget_kw("DX")[0]),
            )
        else:
            dil["arat_array"][dig["actind"]] = np.divide(
                np.array(dig["init"].iget_kw("DZ")[0]),
                (
                    np.array(dig["init"].iget_kw("DX")[0]) ** 2
                    + np.array(dig["init"].iget_kw("DY")[0]) ** 2
                )
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
    if dig["use"] == "opm":
        dil["co2mb_array"][dig["actind"]] = np.array(dig["unrst"]["RES_GAS", t_n + 1])
        if dig["unrst"].count("RES_WAT", t_n + 1):
            dil["h2omb_array"][dig["actind"]] = np.array(
                dig["unrst"]["RES_WAT", t_n + 1]
            )
        else:
            dil["h2omb_array"][dig["actind"]] = np.array(
                dig["unrst"]["RES_OIL", t_n + 1]
            )
    else:
        dil["co2mb_array"][dig["actind"]] = np.array(dig["unrst"]["RES_GAS"][t_n])
        if dig["unrst"].has_kw("RES_WAT"):
            dil["h2omb_array"][dig["actind"]] = np.array(dig["unrst"]["RES_WAT"][t_n])
        else:
            dil["h2omb_array"][dig["actind"]] = np.array(dig["unrst"]["RES_OIL"][t_n])
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
        if dig["case"] == "spe11a":
            dil[f"{name}_array"] = np.empty(dig["nocellst"]) * np.nan
    dil["tco2_array"] = np.zeros(dig["nocellst"])
    dil["tco2_refg"] = np.zeros(dig["nocellsr"])
    dil["tco2_refg"][dig["actindr"]] = np.nan
    if dig["use"] == "opm":
        sgas = abs(np.array(dig["unrst"]["SGAS", t_n]))
        rhog = np.array(dig["unrst"]["GAS_DEN", t_n])
        pres = np.array(dig["unrst"]["PRESSURE", t_n])
        if dig["unrst"].count("PCGW", t_n):
            pres -= np.array(dig["unrst"]["PCGW", t_n])
        rhow = np.array(dig["unrst"][f"{dig['watDen'].upper()}", t_n])
        rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}", t_n])
        if dig["unrst"].count(f"{dig['r_v'].upper()}", t_n):
            rvv = np.array(dig["unrst"][f"{dig['r_v'].upper()}", t_n])
        else:
            rvv = 0.0 * rss
        if dig["case"] != "spe11a":
            dil["temp_array"][dig["actind"]] = np.array(dig["unrst"]["TEMP", t_n])
    else:
        sgas = abs(np.array(dig["unrst"]["SGAS"][t_n]))
        rhog = np.array(dig["unrst"]["GAS_DEN"][t_n])
        pres = np.array(dig["unrst"]["PRESSURE"][t_n])
        if dig["unrst"].has_kw("PCGW"):
            pres -= np.array(dig["unrst"]["PCGW"][t_n])
        rhow = np.array(dig["unrst"][f"{dig['watDen'].upper()}"][t_n])
        rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}"][t_n])
        if dig["unrst"].has_kw(f"{dig['r_v'].upper()}"):
            rvv = np.array(dig["unrst"][f"{dig['r_v'].upper()}"][t_n])
        else:
            rvv = 0.0 * rss
        if dig["case"] != "spe11a":
            dil["temp_array"][dig["actind"]] = np.array(dig["unrst"]["TEMP"][t_n])
    x_l_co2 = np.divide(rss, rss + WAT_DEN_REF / GAS_DEN_REF)
    x_g_h2o = np.divide(rvv, rvv + GAS_DEN_REF / WAT_DEN_REF)
    co2_g = (1 - x_g_h2o) * sgas * rhog * dig["porva"]
    co2_d = x_l_co2 * (1 - sgas) * rhow * dig["porva"]
    dil["pressure_array"][dig["actind"]] = 1e5 * pres
    dil["sgas_array"][dig["actind"]] = sgas * (sgas > SGAS_THR)
    dil["gden_array"][dig["actind"]] = rhog * (sgas > SGAS_THR)
    dil["wden_array"][dig["actind"]] = rhow
    dil["xco2_array"][dig["actind"]] = x_l_co2
    dil["xh20_array"][dig["actind"]] = x_g_h2o * (sgas > SGAS_THR)
    dil["tco2_array"][dig["actind"]] = co2_d + co2_g


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
    for zcord in dil["refzcent"]:
        idxy = 0
        for ycord in dil["refycent"]:
            for xcord in dil["refxcent"]:
                idc = -dig["nxyz"][0] * dig["nxyz"][1] * (dig["nxyz"][2] - idz) + idxy
                if np.isnan(dil["tco2_refg"][idc]):
                    co2 = "n/a"
                else:
                    co2 = f"{dil['tco2_refg'][idc] :.3e}"
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
                            + f"{dil['xco2_refg'][idc] :.3e}, "
                            + f"{dil['xh20_refg'][idc] :.3e}, "
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
                            + f"{dil['xco2_refg'][idc] :.3e}, "
                            + f"{dil['xh20_refg'][idc] :.3e}, "
                            + f"{dil['gden_refg'][idc] :.3e}, "
                            + f"{dil['wden_refg'][idc] :.3e}, "
                            + f"{co2}, "
                            + f"{dil['temp_refg'][idc] :.3e}"
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
                            + f"{dil['xco2_refg'][idc] :.3e}, "
                            + f"{dil['xh20_refg'][idc] :.3e}, "
                            + f"{dil['gden_refg'][idc] :.3e}, "
                            + f"{dil['wden_refg'][idc] :.3e}, "
                            + f"{co2}, "
                            + f"{dil['temp_refg'][idc] :.3e}"
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
