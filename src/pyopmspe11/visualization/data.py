# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT
# pylint: disable=C0302, R0912, R0914

""""
Script to write the benchmark data
"""

import os
import argparse
import csv
from io import StringIO
from shapely.geometry import Polygon
from rtree import index
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

try:
    from opm.io.ecl import EclFile as OpmFile
    from opm.io.ecl import EGrid as OpmGrid
    from opm.io.ecl import ERst as OpmRestart
    from opm.io.ecl import ESmry as OpmSummary
except ImportError:
    print("The Python package opm was not found, using resdata")
try:
    from resdata.grid import Grid
    from resdata.resfile import ResdataFile
    from resdata.summary import Summary
except ImportError:
    print("The resdata Python package was not found, using opm")

GAS_DEN_REF = 1.86843
WAT_DEN_REF = 998.108
SECONDS_IN_YEAR = 31536000
KMOL_TO_KG = 1e3 * 0.044


def main():
    """Postprocessing"""
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
        default="25",
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
        "'dense_sparse', or 'all'",
    )
    parser.add_argument(
        "-u",
        "--use",
        default="resdata",
        help="Using the 'resdata' or python package (resdata by default).",
    )
    cmdargs = vars(parser.parse_known_args()[0])
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
    dig = read_times(dig)
    if dig["use"] == "opm":
        dig = read_opm(dig)
    else:
        dig = read_resdata(dig)
    if dig["mode"] in ["performance", "all", "dense_performance", "performance_sparse"]:
        performance(dig)
    if dig["mode"] in ["all", "sparse", "dense_sparse", "performance_sparse"]:
        sparse_data(dig)
    if dig["mode"] in [
        "all",
        "performance-spatial",
        "dense",
        "dense_performance",
        "dense_sparse",
        "dense_performance-spatial",
    ]:
        if isinstance(dig["dense_t"], float):
            dig["dense_t"] = [
                i * dig["dense_t"]
                for i in range(int(np.floor((dig["times"][-1]) / dig["dense_t"])) + 1)
            ]
        dense_data(dig)


def read_times(dig):
    """Get the time for injection and restart number"""
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
    return dig


def read_resdata(dig):
    """Using resdata"""
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
    return dig


def read_opm(dig):
    """Using opm"""
    dig["unrst"] = OpmRestart(f"{dig['sim']}.UNRST")
    dig["init"] = OpmFile(f"{dig['sim']}.INIT")
    dig["egrid"] = OpmGrid(f"{dig['sim']}.EGRID")
    dig["smspec"] = OpmSummary(f"{dig['sim']}.SMSPEC")
    dig = opm_files(dig)
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
    return dig


def opm_files(dig):
    """Extract the data from the opm output files"""
    dig["norst"] = len(dig["unrst"].report_steps)
    if dig["unrst"].count("WAT_DEN", 0):
        dig["watDen"], dig["r_s"], dig["r_v"] = "wat_den", "rsw", "rvw"
        dig["bpr"] = "BWPR"
    else:
        dig["watDen"], dig["r_s"], dig["r_v"] = "oil_den", "rs", "rv"
        dig["bpr"] = "BPR"
    dig["porv"] = np.array(dig["init"]["PORV"])
    return dig


def performance(dig):
    """Write the performance within the benchmark format"""
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
    # dil["runtimes"] = np.array(
    #     [sum(infostep[i] for i in [2, 3, 4, 5, 6]) for infostep in dil["infosteps"]]
    # )
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
    """Write the performance data"""
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
    for j, time in enumerate(dil["times_data"]):
        ind = j == dil["map_info"]
        itd = j == dil["map_sum"]
        dil["tstep"] = np.sum(dil["tsteps"][ind])
        if dil["tstep"] > 0:
            dil["tstep"] /= np.sum(ind)
        dil["text"].append(
            f"{time:.3e}, "
            + f"{dil['tstep']:.3e}, "
            + f"{np.sum(dil['fsteps'][ind]):.3e}, "
            + f"{interp_fgip(time):.3e}, "
            + f"{dig['dof'] * dig['nocellsa']:.3e}, "
            + f"{np.sum(dil['nliters'][ind]):.3e}, "
            + f"{np.sum(dil['nress'][ind]):.3e}, "
            + f"{np.sum(dil['liniters'][ind]):.3e}, "
            + f"{np.sum(tcpu[itd]):.3e}, "
            + f"{np.sum(dil['tlinsols'][ind]):.3e}"
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
    """Use the summary arrays for the sparse data interpolation"""
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
    if "RGKDI:1" in dig["smspec"].keys():
        dig["X"] = "K"
    else:
        dig["X"] = "C"
    sort = sorted(range(len(i_jk)), key=i_jk.__getitem__)
    if dig["use"] == "opm":
        pop1 = dig["unrst"]["PRESSURE", 0][dil["fipnum"].index(8)]
        pop2 = dig["unrst"]["PRESSURE", 0][dil["fipnum"].index(9)]
        if dig["unrst"].count("PCGW", 0):
            pop1 -= dig["unrst"]["PCGW", 0][dil["fipnum"].index(8)]
            pop2 -= dig["unrst"]["PCGW", 0][dil["fipnum"].index(9)]
        dil["pop1"] = [pop1 * 1.0e5] + list(dig["smspec"][names[sort[0]]] * 1.0e5)  # Pa
        dil["pop2"] = [pop2 * 1.0e5] + list(dig["smspec"][names[sort[1]]] * 1.0e5)  # Pa
        for i in [2, 4, 5, 8]:
            dil["moba"] += dig["smspec"][f"RG{dig['X']}DM:{i}"] * KMOL_TO_KG
            dil["imma"] += dig["smspec"][f"RG{dig['X']}DI:{i}"] * KMOL_TO_KG
            dil["dissa"] += dig["smspec"][f"RWCD:{i}"] * KMOL_TO_KG
        for i in [5, 8]:
            dil["seala"] += (
                dig["smspec"][f"RWCD:{i}"]
                + dig["smspec"][f"RG{dig['X']}DM:{i}"]
                + dig["smspec"][f"RG{dig['X']}DI:{i}"]
            ) * KMOL_TO_KG
        for i in [3, 6]:
            dil["mobb"] += dig["smspec"][f"RG{dig['X']}DM:{i}"] * KMOL_TO_KG
            dil["immb"] += dig["smspec"][f"RG{dig['X']}DI:{i}"] * KMOL_TO_KG
            dil["dissb"] += dig["smspec"][f"RWCD:{i}"] * KMOL_TO_KG
        for key in ["RWCD:6", f"RG{dig['X']}DM:6", f"RG{dig['X']}DI:6"]:
            dil["sealb"] += dig["smspec"][key] * KMOL_TO_KG
        dil["sealt"] = dil["seala"] + dil["sealb"]
        for name in ["RWCD", f"RG{dig['X']}DM", f"RG{dig['X']}DI"]:
            dil["sealt"] += (
                dig["smspec"][f"{name}:7"] + dig["smspec"][f"{name}:9"]
            ) * KMOL_TO_KG
        if dig["case"] != "spe11a":
            sealbound = (
                dig["smspec"]["RWCD:10"]
                + dig["smspec"][f"RG{dig['X']}DM:10"]
                + dig["smspec"][f"RG{dig['X']}DI:10"]
            ) * KMOL_TO_KG
            dil["sealt"] += sealbound
            dil["boundtot"] = (
                sealbound
                + (
                    dig["smspec"]["RWCD:11"]
                    + dig["smspec"][f"RG{dig['X']}DM:11"]
                    + dig["smspec"][f"RG{dig['X']}DI:11"]
                )
                * KMOL_TO_KG
            )
    else:
        dil = resdata_summary(dig, dil, names, sort)
    return dil


def resdata_summary(dig, dil, names, sort):
    """Using resdata"""
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
    for i in [2, 4, 5, 8]:
        dil["moba"] += dig["smspec"][f"RG{dig['X']}DM:{i}"].values * KMOL_TO_KG
        dil["imma"] += dig["smspec"][f"RG{dig['X']}DI:{i}"].values * KMOL_TO_KG
        dil["dissa"] += dig["smspec"][f"RWCD:{i}"].values * KMOL_TO_KG
    for i in [5, 8]:
        dil["seala"] += (
            dig["smspec"][f"RWCD:{i}"].values
            + dig["smspec"][f"RG{dig['X']}DM:{i}"].values
            + dig["smspec"][f"RG{dig['X']}DI:{i}"].values
        ) * KMOL_TO_KG
    for i in [3, 6]:
        dil["mobb"] += dig["smspec"][f"RG{dig['X']}DM:{i}"].values * KMOL_TO_KG
        dil["immb"] += dig["smspec"][f"RG{dig['X']}DI:{i}"].values * KMOL_TO_KG
        dil["dissb"] += dig["smspec"][f"RWCD:{i}"].values * KMOL_TO_KG
    for key in ["RWCD:6", f"RG{dig['X']}DM:6", f"RG{dig['X']}DI:6"]:
        dil["sealb"] += dig["smspec"][key].values * KMOL_TO_KG
    dil["sealt"] = dil["seala"] + dil["sealb"]
    for name in ["RWCD", f"RG{dig['X']}DM", f"RG{dig['X']}DI"]:
        dil["sealt"] += (
            dig["smspec"][f"{name}:7"].values + dig["smspec"][f"{name}:9"].values
        ) * KMOL_TO_KG
    if dig["case"] != "spe11a":
        sealbound = (
            dig["smspec"]["RWCD:10"].values
            + dig["smspec"][f"RG{dig['X']}DM:10"].values
            + dig["smspec"][f"RG{dig['X']}DI:10"].values
        ) * KMOL_TO_KG
        dil["sealt"] += sealbound
        dil["boundtot"] = (
            sealbound
            + (
                dig["smspec"]["RWCD:11"].values
                + dig["smspec"][f"RG{dig['X']}DM:11"].values
                + dig["smspec"][f"RG{dig['X']}DI:11"].values
            )
            * KMOL_TO_KG
        )
    return dil


def overlapping_c_and_facie1_contribution(dig, dil):
    """Add the corresponding fipnum 12 contribution"""
    if dig["use"] == "opm":
        dil["moba"] += dig["smspec"][f"RG{dig['X']}DM:12"] * KMOL_TO_KG
        dil["imma"] += dig["smspec"][f"RG{dig['X']}DI:12"] * KMOL_TO_KG
        dil["dissa"] += dig["smspec"]["RWCD:12"] * KMOL_TO_KG
        dil["seala"] += (
            dig["smspec"]["RWCD:12"]
            + dig["smspec"][f"RG{dig['X']}DM:12"]
            + dig["smspec"][f"RG{dig['X']}DI:12"]
        ) * KMOL_TO_KG
        dil["sealt"] += (
            dig["smspec"]["RWCD:12"]
            + dig["smspec"][f"RG{dig['X']}DM:12"]
            + dig["smspec"][f"RG{dig['X']}DI:12"]
        ) * KMOL_TO_KG
    else:
        dil["moba"] += dig["smspec"][f"RG{dig['X']}DM:12"].values * KMOL_TO_KG
        dil["imma"] += dig["smspec"][f"RG{dig['X']}DI:12"].values * KMOL_TO_KG
        dil["dissa"] += dig["smspec"]["RWCD:12"].values * KMOL_TO_KG
        dil["seala"] += (
            dig["smspec"]["RWCD:12"].values
            + dig["smspec"][f"RG{dig['X']}DM:12"].values
            + dig["smspec"][f"RG{dig['X']}DI:12"].values
        ) * KMOL_TO_KG
        dil["sealt"] += (
            dig["smspec"]["RWCD:12"].values
            + dig["smspec"][f"RG{dig['X']}DM:12"].values
            + dig["smspec"][f"RG{dig['X']}DI:12"].values
        ) * KMOL_TO_KG
    return dil


def sparse_data(dig):
    """Compute the quantities in boxes A, B, and C"""
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
    dil = create_from_summary(dig, dil)
    if dig["case"] == "spe11c" and max(dil["fipnum"]) == 12:
        dil = overlapping_c_and_facie1_contribution(dig, dil)
    # Using the restart data until implemented in OPM Flow summary
    dil = compute_m_c(dig, dil)
    write_sparse_data(dig, dil)


def compute_m_c(dig, dil):
    """Normalized total variation of the concentration field within Box C"""
    dil["boxc"] = np.array([fip in (4, 12) for fip in dil["fipnum"]])
    dil["boxc_x"] = np.roll(dil["boxc"], 1)
    dil["boxc_y"] = np.roll(dil["boxc"], -dig["gxyz"][0])
    dil["boxc_z"] = np.roll(dil["boxc"], -dig["gxyz"][0] * dig["gxyz"][1])
    dil = max_xcw(dig, dil)
    for t_n in range(dig["no_skip_rst"] + 1, dig["norst"]):
        if dig["use"] == "opm":
            sgas = abs(np.array(dig["unrst"]["SGAS", t_n]))
            rhow = np.array(dig["unrst"][f"{dig['watDen'].upper()}", t_n])
            rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}", t_n])
        else:
            sgas = abs(np.array(dig["unrst"]["SGAS"][t_n]))
            rhow = np.array(dig["unrst"][f"{dig['watDen'].upper()}"][t_n])
            rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}"][t_n])
        co2_d = rss * rhow * (1.0 - sgas) * dig["porva"] * GAS_DEN_REF / WAT_DEN_REF
        h2o_l = (1 - sgas) * rhow * dig["porva"]
        mliq = co2_d + h2o_l
        xco2 = 0.0 * co2_d
        inds = mliq > 0.0
        xco2[inds] = np.divide(co2_d[inds], mliq[inds])
        dil["xcw"] = xco2
        if dil["xcw_max"] != 0:
            dil["xcw"] /= dil["xcw_max"]
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
    return dil


def write_sparse_data(dig, dil):
    """Routine to write the sparse data"""
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
                f"{dil['moba'][j]:.3e}, {dil['imma'][j]:.3e}, {dil['dissb'][j]:.3e}, "
                f"{dil['seala'][j]:.3e}, {dil['mobb'][j]:.3e}, {dil['immb'][j]:.3e}, "
                f"{dil['dissb'][j]:.3e}, {dil['sealb'][j]:.3e}, {dil['m_c'][j]:.3e}, "
                f"{dil['sealt'][j]:.3e}"
            )
    else:
        text[-1] += ", boundTot [kg]"
        for j, time in enumerate(dil["times_data"]):
            text.append(
                f"{time:.4e}, {dil['pop1'][j]:.3e}, {dil['pop2'][j]:.3e}, "
                f"{dil['moba'][j]:.3e}, {dil['imma'][j]:.3e}, {dil['dissb'][j]:.3e}, "
                f"{dil['seala'][j]:.3e}, {dil['mobb'][j]:.3e}, {dil['immb'][j]:.3e}, "
                f"{dil['dissb'][j]:.3e}, {dil['sealb'][j]:.3e}, {dil['m_c'][j]:.3e}, "
                f"{dil['sealt'][j]:.3e}, {dil['boundtot'][j]:.3e}"
            )
    with open(
        f"{dig['where']}/{dig['case']}_time_series.csv", "w", encoding="utf8"
    ) as file:
        file.write("\n".join(text))


def max_xcw(dig, dil):
    """Get the maximum CO2 mass fraction in the liquid phase"""
    dil["xcw_max"] = 0
    for t_n in range(dig["no_skip_rst"], dig["norst"]):
        if dig["use"] == "opm":
            sgas = abs(np.array(dig["unrst"]["SGAS", t_n]))
            rhow = np.array(dig["unrst"][f"{dig['watDen'].upper()}", t_n])
            rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}", t_n])
        else:
            sgas = abs(np.array(dig["unrst"]["SGAS"][t_n]))
            rhow = np.array(dig["unrst"][f"{dig['watDen'].upper()}"][t_n])
            rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}"][t_n])
        co2_d = rss * rhow * (1.0 - sgas) * dig["porva"] * GAS_DEN_REF / WAT_DEN_REF
        h2o_l = (1 - sgas) * rhow * dig["porva"]
        mliq = co2_d + h2o_l
        xcw = 0.0 * co2_d
        inds = mliq > 0.0
        xcw[inds] = np.divide(co2_d[inds], mliq[inds])
        xcw_max = np.max(xcw[dil["boxc"]])
        dil["xcw_max"] = max(xcw_max, dil["xcw_max"])
    return dil


def get_corners(dig, dil):
    """Corners from the simulation grid"""
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
    return dil


def dense_data(dig):
    """
    Write the quantities with the benchmark format for the dense data.
    Still plenty of room to improve here the performance and memory usage.
    """
    dil = {"rstno": []}
    for time in dig["dense_t"]:
        dil["rstno"].append(dig["times"].index(time))
    dil = get_corners(dig, dil)
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
        dil = handle_inactive_mapping(dig, dil)
    if dig["case"] == "spe11c":
        dil = handle_yaxis_mapping_intensive(dig, dil)
        dil = handle_yaxis_mapping_extensive(dig, dil)
    if dig["mode"] == "all" or dig["mode"][:5] == "dense":
        names = ["pressure", "sgas", "xco2", "xh20", "gden", "wden", "tco2"]
        if dig["case"] != "spe11a":
            names = ["temp"] + names
        for i, rst in enumerate(dil["rstno"]):
            print(f"Processing dense data {i+1} out of {dil['nrstno']}")
            t_n = rst + dig["no_skip_rst"]
            dil = generate_arrays(dig, dil, names, t_n)
            dil = map_to_report_grid(dig, dil, names)
            write_dense_data(dig, dil, i)
    if dig["mode"] in ["all", "performance-spatial", "dense_performance-spatial"]:
        handle_performance_spatial(dig, dil)


def handle_yaxis_mapping_extensive(dig, dil):
    """Extend the indices accounting for the y direction (extensive quantities)"""
    simycent = [0.0] * dig["gxyz"][1]
    with open(f"{dig['path']}/deck/ycenters.txt", "r", encoding="utf8") as file:
        for j, row in enumerate(csv.reader(file)):
            simycent[j] = float(row[0])
    indy = np.array(
        [pd.Series(np.abs(dil["refycent"] - y_c)).argmin() for y_c in simycent]
    )
    weights = [1.0 / (np.sum(indy == val)) for val in indy]
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
                    col[1] * weights[0],
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
                    [col[0] + (iy * dig["nxyz"][0]), col[1] * weights[j + 1]]
                    for col in row
                ]
                for row in maps
            ]
    dil["cell_ind"] = wei_inds
    return dil


def handle_yaxis_mapping_intensive(dig, dil):
    """Extend the indices accounting for the y direction"""
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
    return dil


def handle_inactive_mapping(dig, dil):
    """Set to inf the inactive grid centers in the reporting grid"""
    for i in dig["actind"]:
        for mask in dil["cell_ind"][i]:
            dig["actindr"].append(mask[0])
    dig["actindr"] = list(dict.fromkeys(dig["actindr"]))
    allc = np.linspace(0, dig["nocellsr"] - 1, dig["nocellsr"], dtype=int)
    dig["actindr"] = np.delete(allc, dig["actindr"])
    return dil


def handle_performance_spatial(dig, dil):
    """Create the performance spatial maps"""
    dil["counter"] = 0.0 * np.ones(dig["nocellsr"])
    dil["pv"] = 0.0 * np.ones(dig["nocellsr"])
    dil["pv"][dig["actindr"]] = 1.0
    dil = static_map_to_report_grid_performance_spatial(dig, dil)
    names = ["co2mn", "h2omn", "co2mb", "h2omb"]
    for i, rst in enumerate(dil["rstno"]):
        print(f"Processing performance spatial {i+1} out of {dil['nrstno']}")
        for name in names:
            dil[f"{name}_array"] = np.zeros(dig["nocellst"])
        t_n = rst + dig["no_skip_rst"]
        if t_n > 0:
            # RESIDUAL not included in the SOLUTION deck section (substract 1)
            dil = generate_arrays_performance_spatial(dig, dil, t_n - 1)
        dil = map_to_report_grid_performance_spatial(
            dig, dil, names, dil["latest_dts"][i]
        )
        write_dense_data_performance_spatial(dig, dil, i)


def static_map_to_report_grid_performance_spatial(dig, dil):
    """Map the no dynamic quantities to the reporting grid"""
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
    return dil


def generate_arrays_performance_spatial(dig, dil, t_n):
    """Numpy arrays for the performance spatial data"""
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
    return dil


def map_to_report_grid_performance_spatial(dig, dil, names, d_t):
    """Map the simulation grid to the reporting grid"""
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
    return dil


def write_dense_data_performance_spatial(dig, dil, i):
    """Generate the cvs"""
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
    """Numpy arrays for the dense data"""
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
    co2_g = sgas * rhog * dig["porva"]
    co2_d = rss * rhow * (1.0 - sgas) * dig["porva"] * GAS_DEN_REF / WAT_DEN_REF
    h2o_l = (1 - sgas) * rhow * dig["porva"]
    dil["pressure_array"][dig["actind"]] = 1e5 * pres
    dil["sgas_array"][dig["actind"]] = sgas
    dil["gden_array"][dig["actind"]] = rhog * (sgas > 0.0)
    dil["wden_array"][dig["actind"]] = rhow
    dil["tco2_array"][dig["actind"]] = co2_d + co2_g
    dil = compute_xco2(dig, dil, co2_d, h2o_l)
    h2o_v = rvv * rhog * sgas * dig["porva"] * WAT_DEN_REF / GAS_DEN_REF
    dil = compute_xh20(dig, dil, h2o_v, co2_g)
    return dil


def compute_xco2(dig, dil, co2_d, h2o_l):
    """Compute the mass fraction of CO2 in liquid"""
    mliq = co2_d + h2o_l
    xco2 = 0.0 * co2_d
    inds = mliq > 0.0
    xco2[inds] = np.divide(co2_d[inds], mliq[inds])
    dil["xco2_array"][dig["actind"]] = xco2
    return dil


def compute_xh20(dig, dil, h2o_v, co2_g):
    """Compute the mass fraction of water in vapour"""
    mgas = h2o_v + co2_g
    xh20 = 0.0 * h2o_v
    inds = mgas > 0.0
    xh20[inds] = np.divide(h2o_v[inds], mgas[inds])
    dil["xh20_array"][dig["actind"]] = xh20
    return dil


def map_to_report_grid(dig, dil, names):
    """Map the simulation grid to the reporting grid"""
    for i in dig["actind"]:
        for mask in dil["cell_ind"][i]:
            dil["tco2_refg"][mask[0]] += dil["tco2_array"][i] * mask[1]
    for i, ind in enumerate(dil["cell_cent"]):
        for name in names[:-1]:
            dil[f"{name}_refg"][i] = dil[f"{name}_array"][int(ind)]
    return dil


def write_dense_data(dig, dil, i):
    """Map the quantities to the cells"""
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
    """Get the right file header"""
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
