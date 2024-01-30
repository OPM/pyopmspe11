# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

""""
Script to write the benchmark data
"""

import os
import argparse
import csv
from io import StringIO
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
    parser = argparse.ArgumentParser(description="Main script to porcess the data")
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
        help="Write only the 'dense', 'performance', 'sparse', 'dense_performance', "
        "'dense_sparse', 'performance_sparse', or 'all'",
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
    dig = read_times(dig)
    if dig["use"] == "opm":
        dig = read_opm(dig)
    else:
        dig = read_resdata(dig)
    if dig["mode"] in ["performance", "all", "dense_performance", "performance_sparse"]:
        performance(dig)
    if dig["mode"] in ["all", "sparse", "dense_sparse", "performance_sparse"]:
        sparse_data(dig)
    if dig["mode"] in ["all", "dense", "dense_performance", "dense_sparse"]:
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
    """Write the performance within the benchmark format SECONDS_IN_YEAR"""
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
    times = max(0, dig["no_skip_rst"] - 1)
    dil["map_info"] = np.array(
        [times + int(np.floor(infotime / dig["sparse_t"])) for infotime in infotimes]
    )
    dil["fsteps"] = np.array(
        [1.0 * (infostep[11] == 0) for infostep in dil["infosteps"]]
    )
    dil["nress"] = np.array([infostep[8] for infostep in dil["infosteps"]])
    dil["tlinsols"] = np.array([infostep[4] for infostep in dil["infosteps"]])
    dil["runtimes"] = np.array(
        [sum(infostep[i] for i in [2, 4, 5, 6]) for infostep in dil["infosteps"]]
    )
    dil["liniters"] = np.array([infostep[10] for infostep in dil["infosteps"]])
    dil["nliters"] = np.array([infostep[9] for infostep in dil["infosteps"]])
    dil["tsteps"] = np.array(
        [86400 * infostep[1] * infostep[11] for infostep in dil["infosteps"]]
    )
    if dig["use"] == "opm":
        fgip = dig["smspec"]["FGIP"]
        times = 86400.0 * dig["smspec"]["TIME"] - dig["time_initial"]
    else:
        fgip = dig["smspec"]["FGIP"].values
        times = 86400.0 * dig["smspec"]["TIME"].values - dig["time_initial"]
    interp_fgip = interp1d(
        times,
        GAS_DEN_REF * fgip,
        fill_value="extrapolate",
    )
    dil["text"] = []
    dil["text"].append(
        "# t [s], tstep [s], fsteps [-], mass [kg], dof [-], nliter [-], "
        + "nres [-], liniter [-], runtime [s], tlinsol [s]"
    )
    if dig["no_skip_rst"] == 0:
        dil["text"].append(
            "0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, "
            + "0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00"
        )
        dil["times_data"] = np.delete(dil["times_data"], 0)
    for j, time in enumerate(dil["times_data"]):
        ind = j == dil["map_info"]
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
            + f"{np.sum(dil['runtimes'][ind]):.3e}, "
            + f"{np.sum(dil['tlinsols'][ind]):.3e}"
        )
    with open(
        f"{dig['where']}/{dig['case']}_performance_time_series.csv",
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
            dil["moba"] += dig["smspec"][f"RGCDM:{i}"] * KMOL_TO_KG
            dil["imma"] += dig["smspec"][f"RGCDI:{i}"] * KMOL_TO_KG
            dil["dissa"] += dig["smspec"][f"RWCD:{i}"] * KMOL_TO_KG
        for i in [5, 8]:
            dil["seala"] += (
                dig["smspec"][f"RWCD:{i}"]
                + dig["smspec"][f"RGCDM:{i}"]
                + dig["smspec"][f"RGCDI:{i}"]
            ) * KMOL_TO_KG
        for i in [3, 6]:
            dil["mobb"] += dig["smspec"][f"RGCDM:{i}"] * KMOL_TO_KG
            dil["immb"] += dig["smspec"][f"RGCDI:{i}"] * KMOL_TO_KG
            dil["dissb"] += dig["smspec"][f"RWCD:{i}"] * KMOL_TO_KG
        for key in ["RWCD:6", "RGCDM:6", "RGCDI:6"]:
            dil["sealb"] += dig["smspec"][key] * KMOL_TO_KG
        dil["sealt"] = dil["seala"] + dil["sealb"]
        for name in ["RWCD", "RGCDM", "RGCDI"]:
            dil["sealt"] += (
                dig["smspec"][f"{name}:7"] + dig["smspec"][f"{name}:9"]
            ) * KMOL_TO_KG
        if dig["case"] != "spe11a":
            sealbound = (
                dig["smspec"]["RWCD:10"]
                + dig["smspec"]["RGCDM:10"]
                + dig["smspec"]["RGCDI:10"]
            ) * KMOL_TO_KG
            dil["sealt"] += sealbound
            dil["boundtot"] = (
                sealbound
                + (
                    dig["smspec"]["RWCD:11"]
                    + dig["smspec"]["RGCDM:11"]
                    + dig["smspec"]["RGCDI:11"]
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
        dil["moba"] += dig["smspec"][f"RGCDM:{i}"].values * KMOL_TO_KG
        dil["imma"] += dig["smspec"][f"RGCDI:{i}"].values * KMOL_TO_KG
        dil["dissa"] += dig["smspec"][f"RWCD:{i}"].values * KMOL_TO_KG
    for i in [5, 8]:
        dil["seala"] += (
            dig["smspec"][f"RWCD:{i}"].values
            + dig["smspec"][f"RGCDM:{i}"].values
            + dig["smspec"][f"RGCDI:{i}"].values
        ) * KMOL_TO_KG
    for i in [3, 6]:
        dil["mobb"] += dig["smspec"][f"RGCDM:{i}"].values * KMOL_TO_KG
        dil["immb"] += dig["smspec"][f"RGCDI:{i}"].values * KMOL_TO_KG
        dil["dissb"] += dig["smspec"][f"RWCD:{i}"].values * KMOL_TO_KG
    for key in ["RWCD:6", "RGCDM:6", "RGCDI:6"]:
        dil["sealb"] += dig["smspec"][key].values * KMOL_TO_KG
    dil["sealt"] = dil["seala"] + dil["sealb"]
    for name in ["RWCD", "RGCDM", "RGCDI"]:
        dil["sealt"] += (
            dig["smspec"][f"{name}:7"].values + dig["smspec"][f"{name}:9"].values
        ) * KMOL_TO_KG
    if dig["case"] != "spe11a":
        sealbound = (
            dig["smspec"]["RWCD:10"].values
            + dig["smspec"]["RGCDM:10"].values
            + dig["smspec"]["RGCDI:10"].values
        ) * KMOL_TO_KG
        dil["sealt"] += sealbound
        dil["boundtot"] = (
            sealbound
            + (
                dig["smspec"]["RWCD:11"].values
                + dig["smspec"]["RGCDM:11"].values
                + dig["smspec"]["RGCDI:11"].values
            )
            * KMOL_TO_KG
        )
    return dil


def overlapping_c_and_facie1_contribution(dig, dil):
    """Add the corresponding fipnum 12 contribution"""
    if dig["use"] == "opm":
        dil["moba"] += dig["smspec"]["RGCDM:12"] * KMOL_TO_KG
        dil["imma"] += dig["smspec"]["RGCDI:12"] * KMOL_TO_KG
        dil["dissa"] += dig["smspec"]["RWCD:12"] * KMOL_TO_KG
        dil["seala"] += (
            dig["smspec"]["RWCD:12"]
            + dig["smspec"]["RGCDM:12"]
            + dig["smspec"]["RGCDI:12"]
        ) * KMOL_TO_KG
        dil["sealt"] += (
            dig["smspec"]["RWCD:12"]
            + dig["smspec"]["RGCDM:12"]
            + dig["smspec"]["RGCDI:12"]
        ) * KMOL_TO_KG
    else:
        dil["moba"] += dig["smspec"]["RGCDM:12"].values * KMOL_TO_KG
        dil["imma"] += dig["smspec"]["RGCDI:12"].values * KMOL_TO_KG
        dil["dissa"] += dig["smspec"]["RWCD:12"].values * KMOL_TO_KG
        dil["seala"] += (
            dig["smspec"]["RWCD:12"].values
            + dig["smspec"]["RGCDM:12"].values
            + dig["smspec"]["RGCDI:12"].values
        ) * KMOL_TO_KG
        dil["sealt"] += (
            dig["smspec"]["RWCD:12"].values
            + dig["smspec"]["RGCDM:12"].values
            + dig["smspec"]["RGCDI:12"].values
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
            sgas = np.array(dig["unrst"]["SGAS", t_n])
            rhow = np.array(dig["unrst"][f"{dig['watDen'].upper()}", t_n])
            rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}", t_n])
        else:
            sgas = np.array(dig["unrst"]["SGAS"][t_n])
            rhow = np.array(dig["unrst"][f"{dig['watDen'].upper()}"][t_n])
            rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}"][t_n])
        co2_d = rss * rhow * (1.0 - sgas) * dig["porva"] * GAS_DEN_REF / WAT_DEN_REF
        h2o_l = (1 - sgas) * rhow * dig["porva"]
        dil["xcw"] = np.divide(co2_d, co2_d + h2o_l)
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
    if dig["case"] != "spe11a":
        text[-1] += ", boundTot [kg]"
    for j, time in enumerate(dil["times_data"]):
        text.append(
            f"{time:.3e},{dil['pop1'][j]:.3e},{dil['pop2'][j]:.3e}"
            f",{dil['moba'][j]:.3e},{dil['imma'][j]:.3e},{dil['dissb'][j]:.3e}"
            f",{dil['seala'][j]:.3e},{dil['mobb'][j]:.3e},{dil['immb'][j]:.3e}"
            f",{dil['dissb'][j]:.3e},{dil['sealb'][j]:.3e},{dil['m_c'][j]:.3e}"
            f",{dil['sealt'][j]:.3e}"
        )
        if dig["case"] != "spe11a":
            text[-1] += f",{dil['boundtot'][j]:.3e}"
    with open(
        f"{dig['where']}/{dig['case']}_time_series.csv", "w", encoding="utf8"
    ) as file:
        file.write("\n".join(text))


def max_xcw(dig, dil):
    """Get the maximum CO2 mass fraction in the liquid phase"""
    dil["xcw_max"] = 0
    for t_n in range(dig["no_skip_rst"], dig["norst"]):
        if dig["use"] == "opm":
            sgas = np.array(dig["unrst"]["SGAS", t_n])
            rhow = np.array(dig["unrst"][f"{dig['watDen'].upper()}", t_n])
            rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}", t_n])
        else:
            sgas = np.array(dig["unrst"]["SGAS"][t_n])
            rhow = np.array(dig["unrst"][f"{dig['watDen'].upper()}"][t_n])
            rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}"][t_n])
        co2_d = rss * rhow * (1.0 - sgas) * dig["porva"] * GAS_DEN_REF / WAT_DEN_REF
        h2o_l = (1 - sgas) * rhow * dig["porva"]
        xcw = np.divide(co2_d, co2_d + h2o_l)
        xcw_max = np.max(xcw[dil["boxc"]])
        dil["xcw_max"] = max(xcw_max, dil["xcw_max"])
    return dil


def get_centers(dig, dil):
    """Centers from the simulation grid"""
    for i in ["x", "y", "z"]:
        dil[f"sim{i}cent"] = [0.0] * dig["nocellst"]
    if dig["use"] == "opm":
        # opm.io.ecl.EGrid.xyz_from_ijk is returning nan values, that's why we do this
        with open(f"{dig['path']}/deck/centers.txt", "r", encoding="utf8") as file:
            for j, row in enumerate(csv.reader(file)):
                dil["simxcent"][j] = float(row[0])
                dil["simycent"][j] = float(row[1])
                dil["simzcent"][j] = dig["dims"][2] - float(row[2])
        dil["satnum"] = list(dig["init"]["SATNUM"])
    else:
        for cell in dig["egrid"].cells():
            (
                dil["simxcent"][cell.global_index],
                dil["simycent"][cell.global_index],
                dil["simzcent"][cell.global_index],
            ) = (
                cell.coordinate[0],
                cell.coordinate[1],
                dig["dims"][2] - cell.coordinate[2],
            )
        dil["satnum"] = list(dig["init"].iget_kw("SATNUM")[0])
    return dil


def dense_data(dig):
    """Write the quantities with the benchmark format"""
    dil = {"rstno": []}
    for time in dig["dense_t"]:
        dil["rstno"].append(dig["times"].index(time))
    dil = get_centers(dig, dil)
    dil["nrstno"] = len(dil["rstno"])
    for i, j, k in zip(["x", "y", "z"], dig["dims"], dig["nxyz"]):
        ind = np.linspace(0, j, k + 1)
        dil[f"ref{i}cent"] = 0.5 * (ind[1:] + ind[:-1])
        dil[f"ref{i}grid"] = np.zeros(dig["nocellsr"])
    ind = 0
    for k in dil["refzcent"]:
        for j in dil["refycent"]:
            for i in dil["refxcent"]:
                dil["refxgrid"][ind] = i
                dil["refygrid"][ind] = j
                dil["refzgrid"][ind] = k
                ind += 1
    ind, dil["cell_ind"] = 0, np.zeros(dig["nocellst"], dtype=int)
    for i, j, k in zip(dil["simxcent"], dil["simycent"], np.flip(dil["simzcent"])):
        dil["cell_ind"][ind] = pd.Series(
            (dil["refxgrid"] - i) ** 2
            + (dil["refygrid"] - j) ** 2
            + (dil["refzgrid"] - k) ** 2
        ).argmin()
        ind += 1
    if max(dil["satnum"]) < 7:
        dil = handle_inactive_mapping(dig, dil)
    names = ["pressure", "sgas", "xco2", "xh20", "gden", "wden", "tco2"]
    if dig["case"] != "spe11a":
        names += ["temp"]
    for i, rst in enumerate(dil["rstno"]):
        print(f"Processing dense data {i+1} out of {dil['nrstno']}")
        t_n = rst + dig["no_skip_rst"]
        dil = generate_arrays(dig, dil, names, t_n)
        dil = map_to_report_grid(dig, dil, names)
        write_dense_data(dig, dil, i)


def generate_arrays(dig, dil, names, t_n):
    """Numpy arrays for the dense data"""
    for name in names:
        dil[f"{name}_array"] = np.zeros(dig["nocellst"])
    if dig["use"] == "opm":
        sgas = np.array(dig["unrst"]["SGAS", t_n])
        rhog = np.array(dig["unrst"]["GAS_DEN", t_n])
        pres = np.array(dig["unrst"]["PRESSURE", t_n])
        if dig["unrst"].count("PCGW", t_n):
            pres -= np.array(dig["unrst"]["PCGW", t_n])
        rhow = np.array(dig["unrst"][f"{dig['watDen'].upper()}", t_n])
        rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}", t_n])
        if dig["case"] != "spe11a":
            rvv = np.array(dig["unrst"][f"{dig['r_v'].upper()}", t_n])
            dil["temp_array"][dig["actind"]] = np.array(dig["unrst"]["TEMP", t_n])
    else:
        sgas = np.array(dig["unrst"]["SGAS"][t_n])
        rhog = np.array(dig["unrst"]["GAS_DEN"][t_n])
        pres = np.array(dig["unrst"]["PRESSURE"][t_n])
        if dig["unrst"].has_kw("PCGW"):
            pres -= np.array(dig["unrst"]["PCGW"][t_n])
        rhow = np.array(dig["unrst"][f"{dig['watDen'].upper()}"][t_n])
        rss = np.array(dig["unrst"][f"{dig['r_s'].upper()}"][t_n])
        if dig["case"] != "spe11a":
            rvv = np.array(dig["unrst"][f"{dig['r_v'].upper()}"][t_n])
            dil["temp_array"][dig["actind"]] = np.array(dig["unrst"]["TEMP"][t_n])
    co2_g = sgas * rhog * dig["porva"]
    co2_d = rss * rhow * (1.0 - sgas) * dig["porva"] * GAS_DEN_REF / WAT_DEN_REF
    h2o_l = (1 - sgas) * rhow * dig["porva"]
    dil["pressure_array"][dig["actind"]] = 1e5 * pres
    dil["sgas_array"][dig["actind"]] = sgas
    dil["gden_array"][dig["actind"]] = rhog * (sgas > 0.0)
    dil["wden_array"][dig["actind"]] = rhow
    dil["xco2_array"][dig["actind"]] = np.divide(co2_d, co2_d + h2o_l)
    dil["tco2_array"][dig["actind"]] = co2_d + co2_g
    if dig["case"] != "spe11a":
        h2o_v = rvv * rhog * sgas * dig["porva"] * WAT_DEN_REF / GAS_DEN_REF
        dil = compute_xh20(dig, dil, h2o_v, co2_g)
    return dil


def compute_xh20(dig, dil, h2o_v, co2_g):
    """Compute the mass fraction of water in vapour"""
    mgas = h2o_v + co2_g
    xh20 = 0.0 * h2o_v
    inds = mgas > 0.0
    xh20[inds] = np.divide(h2o_v[inds], mgas[inds])
    dil["xh20_array"][dig["actind"]] = xh20
    return dil


def handle_inactive_mapping(dig, dil):
    """Set to inf the inactive grid centers in the reporting grid"""
    var_array = np.empty(dig["nocellst"]) * np.nan
    var_array[dig["actind"]] = 0.0
    for i in np.unique(dil["cell_ind"]):
        inds = i == dil["cell_ind"]
        if np.isnan(np.sum(var_array[inds])):
            dil["refxgrid"][i] = np.inf
            dil["refygrid"][i] = np.inf
            dil["refzgrid"][i] = np.inf
    ind, dil["cell_ind"] = 0, np.zeros(dig["nocellst"], dtype=int)
    for i, j, k in zip(dil["simxcent"], dil["simycent"], np.flip(dil["simzcent"])):
        dil["cell_ind"][ind] = pd.Series(
            (dil["refxgrid"] - i) ** 2
            + (dil["refygrid"] - j) ** 2
            + (dil["refzgrid"] - k) ** 2
        ).argmin()
        ind += 1
    return dil


def map_to_report_grid(dig, dil, names):
    """Map the simulation grid to the reporting grid"""
    for name in names:
        dil[f"{name}_refg"] = np.empty(dig["nocellsr"]) * np.nan
        for i in np.unique(dil["cell_ind"]):
            inds = i == dil["cell_ind"]
            p_v = np.sum(dig["porv"][inds])
            if p_v > 0:
                if name == "tco2":
                    dil[f"{name}_refg"][i] = np.sum(dil[f"{name}_array"][inds])
                else:
                    dil[f"{name}_refg"][i] = (
                        np.sum(dil[f"{name}_array"][inds] * dig["porv"][inds]) / p_v
                    )
    return dil


def write_dense_data(dig, dil, i):
    """Map the quantities to the cells"""
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
    idz = 0
    for zcord in dil["refzcent"]:
        idxy = 0
        for ycord in dil["refycent"]:
            for xcord in dil["refxcent"]:
                idc = (
                    dig["nxyz"][0] * dig["nxyz"][1] * (dig["nxyz"][2] - idz - 1) + idxy
                )
                if dig["case"] == "spe11a":
                    text.append(
                        f"{xcord:.3e}, {zcord:.3e}, {dil['pressure_refg'][idc] :.3e}, "
                        + f"{dil['sgas_refg'][idc] :.3e}, {dil['xco2_refg'][idc] :.3e}, "
                        + f"{dil['xh20_refg'][idc] :.3e}, {dil['gden_refg'][idc] :.3e}, "
                        + f"{dil['wden_refg'][idc] :.3e}, {dil['tco2_refg'][idc] :.3e}"
                    )
                elif dig["case"] == "spe11b":
                    text.append(
                        f"{xcord:.3e}, {zcord:.3e}, {dil['pressure_refg'][idc] :.3e}, "
                        + f"{dil['sgas_refg'][idc] :.3e}, {dil['xco2_refg'][idc] :.3e}, "
                        + f"{dil['xh20_refg'][idc] :.3e}, {dil['gden_refg'][idc] :.3e}, "
                        + f"{dil['wden_refg'][idc] :.3e}, {dil['tco2_refg'][idc] :.3e}, "
                        + f"{dil['temp_refg'][idc] :.3e}"
                    )
                else:
                    text.append(
                        f"{xcord:.3e}, {ycord:.3e}, {zcord:.3e}, "
                        + f"{dil['pressure_refg'][idc] :.3e}, {dil['sgas_refg'][idc] :.3e}, "
                        + f"{dil['xco2_refg'][idc] :.3e}, {dil['xh20_refg'][idc] :.3e}, "
                        + f"{dil['gden_refg'][idc] :.3e}, {dil['wden_refg'][idc] :.3e}, "
                        + f"{dil['tco2_refg'][idc] :.3e}, {dil['temp_refg'][idc] :.3e}"
                    )
                idxy += 1
        idz += 1
    with open(
        f"{dig['where']}/{dig['case']}_spatial_map_{name_t}.csv", "w", encoding="utf8"
    ) as file:
        file.write("\n".join(text))


if __name__ == "__main__":
    main()
