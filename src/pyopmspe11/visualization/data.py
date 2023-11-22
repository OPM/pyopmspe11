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
from resdata.grid import Grid
from resdata.resfile import ResdataFile
from resdata.summary import Summary

GAS_DEN_REF = 1.86843
WAT_DEN_REF = 998.108
SECONDS_IN_YEAR = 31536000


def main():
    """Postprocessing"""
    parser = argparse.ArgumentParser(description="Main script to plot the results")
    parser.add_argument(
        "-p",
        "--path",
        default="output",
        help="The name of the output folder.",
    )
    parser.add_argument(
        "-d",
        "--deck",
        default="spe11b",
        help="The simulated case.",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        default="10,10,5",
        help="Number of x, y, and z elements to write the data ('10,10,5' by default).",
    )
    parser.add_argument(
        "-t",
        "--time",
        default="25",
        help="Time interval for the spatial maps (spe11a [h]; spe11b/c [y]) ('24' by default).",
    )
    parser.add_argument(
        "-g",
        "--generate",
        default="sparse",
        help="Write only the 'dense', 'performance', 'sparse', 'dense_performance', "
        "'dense_sparse', 'performance_sparse', or 'all'",
    )
    cmdargs = vars(parser.parse_known_args()[0])
    dic = {"path": cmdargs["path"].strip()}
    dic["case"] = cmdargs["deck"].strip()
    dic["mode"] = cmdargs["generate"].strip()
    dic["exe"] = os.getcwd()
    dic["where"] = f"{dic['exe']}/{dic['path']}/data"
    dic["nxyz"] = np.genfromtxt(
        StringIO(cmdargs["resolution"]), delimiter=",", dtype=int
    )
    if dic["case"] == "spe11a":
        dic["spatial_t"] = float(cmdargs["time"].strip()) * 3600
        dic["dims"] = [2.8, 1.0, 1.2]
        dic["dof"] = 2
        dic["t0_rst"] = 0
        dic["t1_rst"] = 0
        dic["nxyz"][1] = 1
    else:
        dic["spatial_t"] = float(cmdargs["time"].strip()) * SECONDS_IN_YEAR
        dic["dims"] = [8400.0, 1.0, 1200.0]
        dic["dof"] = 3
        dic["t0_rst"] = 1
        dic["t1_rst"] = 2
        dic["nxyz"][1] = 1
    if dic["case"] == "spe11c":
        dic["dims"][1] = 5000.0
    dic["nocellsr"] = dic["nxyz"][0] * dic["nxyz"][1] * dic["nxyz"][2]
    case = "./" + dic["path"] + "/flow/" + f"{dic['path'].upper()}"
    for name in ["unrst", "init"]:
        dic[f"{name}"] = ResdataFile(case + f".{name.upper()}")
    dic["egrid"], dic["smspec"] = Grid(case + ".EGRID"), Summary(case + ".SMSPEC")
    names = ["sgas", "rs", "pressure", "gas_den", "oil_den", "gaskr"]
    if dic["case"] != "spe11a":
        names += ["rv", "temp"]
    for name in names:
        dic[f"{name}"] = list(dic["unrst"].iget_kw(name.upper()))
    for name in ["porv", "satnum", "fipnum", "dx", "dy", "dz"]:
        dic[f"{name}"] = list(dic["init"].iget_kw(name.upper())[0])
    dic["actnum"] = list(dic["egrid"].export_actnum())
    dic["actind"] = list(i for i, act in enumerate(dic["actnum"]) if act == 1)
    dic["porva"] = [porv for (porv, act) in zip(dic["porv"], dic["actnum"]) if act == 1]
    dic["nocellst"], dic["nocellsa"] = len(dic["actnum"]), sum(dic["actnum"])
    dic["norst"] = dic["unrst"].num_report_steps()
    if dic["mode"] in ["performance", "all", "dense_performance", "performance_sparse"]:
        performance(dic)
    if dic["mode"] in ["all", "dense", "dense_performance", "dense_sparse"]:
        dense_data(dic)
    if dic["mode"] in ["all", "sparse", "dense_sparse", "performance_sparse"]:
        sparse_data(dic)


def performance(dic):
    """Write the performance within the benchmark format"""
    dic["d_t"] = (
        dic["smspec"].end_time - dic["unrst"].dates[dic["t0_rst"]]
    ).total_seconds() / (len(dic["unrst"].dates) - 1 - dic["t0_rst"])
    dic["times"] = [dic["d_t"] * j for j in range(1, dic["norst"] - 1 - dic["t0_rst"])]
    dic["t_0"] = (
        dic["unrst"].dates[dic["t0_rst"]] - dic["unrst"].dates[0]
    ).total_seconds()
    dic["infosteps"] = []
    with open(
        f"{dic['path']}/flow/{dic['path'].upper()}.INFOSTEP", "r", encoding="utf8"
    ) as file:
        for j, row in enumerate(csv.reader(file)):
            if j > 0:
                if float((row[0].strip()).split()[0]) >= (dic["t_0"] / 86400.0):
                    dic["infosteps"].append(
                        [float(column) for column in (row[0].strip()).split()]
                    )
    infotimes = [infostep[0] * 86400 - dic["t_0"] for infostep in dic["infosteps"]]
    dic["map_info"] = [round(infotime / dic["d_t"]) for infotime in infotimes]
    fsteps = [infostep[11] for infostep in dic["infosteps"]]
    nress = [infostep[8] for infostep in dic["infosteps"]]
    tlinsols = [infostep[4] for infostep in dic["infosteps"]]
    runtimes = [sum(infostep[2:7]) for infostep in dic["infosteps"]]
    rsteps = [
        rstep for rstep in dic["smspec"].get_report_step() if rstep > dic["t0_rst"]
    ]
    fgip = dic["smspec"]["FGIP"]
    tsteps = dic["smspec"]["TIMESTEP"]
    nliters = dic["smspec"]["NEWTON"]
    liniters = dic["smspec"]["MLINEARS"]
    dic["map_rsteps"] = [
        sum((1 if (rstep < i) else 0 for rstep in dic["smspec"].get_report_step()))
        for i in range(dic["t1_rst"] + 1, dic["norst"] + 1)
    ]
    dic["text"] = []
    dic["text"].append(
        "# t [s], tstep [s], fsteps [-], mass [kg], dof [-], nliter [-], "
        + "nres [-], liniter [-], runtime [s], tlinsol [s]"
    )
    if dic["t0_rst"] == 0:
        dic["text"].append(
            "0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, "
            + "0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00"
        )
    else:
        dic["times"].insert(0, 0.0)
    for j, time in enumerate(dic["times"]):
        dic["tstep"] = sum(
            (
                tsteps[i].value
                for i, rstep in enumerate(rsteps)
                if rstep == (j + 1 + dic["t0_rst"])
            )
        ) / rsteps.count(j + 1 + dic["t0_rst"])
        dic["nliters"] = sum(
            (
                nliters[i].value
                for i, rstep in enumerate(rsteps)
                if rstep == (j + 1 + dic["t0_rst"])
            )
        )
        dic["liniters"] = sum(
            (
                liniters[i].value
                for i, rstep in enumerate(rsteps)
                if rstep == (j + 1 + dic["t0_rst"])
            )
        )
        dic["text"].append(
            f"{time:.3e}, "
            + f"{dic['tstep']:.3e}, "
            + f"{sum((fsteps[i] for i in dic['map_info'] if i==j)):.3e}, "
            + f"{GAS_DEN_REF*fgip[dic['map_rsteps'][j]-1].value:.3e}, "
            + f"{dic['dof'] * dic['nocellsa']:.3e}, "
            + f"{dic['nliters']:.3e}, "
            + f"{sum((nress[i] for i in dic['map_info'] if i==j)):.3e}, "
            + f"{dic['liniters']:.3e}, "
            + f"{sum((runtimes[i] for i in dic['map_info'] if i==j)):.3e}, "
            + f"{sum((tlinsols[i] for i in dic['map_info'] if i==j)):.3e}"
        )
    with open(
        f"{dic['where']}/{dic['case']}_performance_time_series.csv",
        "w",
        encoding="utf8",
    ) as file:
        file.write("\n".join(dic["text"]))


def sparse_data(dic):
    """Compute the quantities in boxes A, B, and C"""
    for ent in [
        "ip1c",
        "ip2c",
        "moba",
        "imma",
        "dissa",
        "seala",
        "mobb",
        "immb",
        "dissb",
        "sealb",
        "m_c",
        "sealtot",
    ]:
        dic[ent] = []
    if dic["case"] != "spe11a":
        dic["boundtot"] = []
    dic["d_t"] = (
        dic["smspec"].end_time - dic["unrst"].dates[dic["t0_rst"]]
    ).total_seconds() / (len(dic["unrst"].dates) - 1 - dic["t0_rst"])
    dic["boxa"] = [i for i, fip in enumerate(dic["fipnum"]) if fip in (2, 4, 5)]
    dic["boxb"] = [i for i, fip in enumerate(dic["fipnum"]) if fip == 3]
    dic["boxc"] = np.array([fip == 4 for fip in dic["fipnum"]])
    dic["boxc_x"] = list(np.roll(dic["boxc"], 1))
    dic["boxc_y"] = list(np.roll(dic["boxc"], -dic["egrid"].nx))
    dic["boxc_z"] = list(np.roll(dic["boxc"], -dic["egrid"].nx * dic["egrid"].ny))
    dic["boxc"] = [i for i, fip in enumerate(dic["fipnum"]) if fip == 4]
    dic["boxc_z"] = [i for i, fip in enumerate(dic["boxc_z"]) if fip == 1]
    dic["boxc_y"] = [i for i, fip in enumerate(dic["boxc_y"]) if fip == 1]
    dic["boxc_x"] = [i for i, fip in enumerate(dic["boxc_x"]) if fip == 1]
    dic["sensor1"] = dic["fipnum"].index(5)
    dic["sensor2"] = dic["fipnum"].index(6)
    dic["facie1"] = [sat == 1 for sat in dic["satnum"]]
    dic["facie1ind"] = [i for i, sat in enumerate(dic["satnum"]) if sat == 1]
    dic["boundariesind"] = [i for i, fip in enumerate(dic["fipnum"]) if fip == 7]
    write_sparse_data(dic)


def write_sparse_data(dic):
    """Map the quantities to the locations"""
    text = [
        "# t [s], p1 [Pa], p2 [Pa], mobA [kg], immA [kg], dissA [kg], sealA [kg], "
        + "<same for B>, MC [m^2], sealTot [kg]"
    ]
    if dic["case"] != "spe11a":
        text[-1] += ", boundTot [kg]"
    for tim, j in enumerate(range(dic["t1_rst"], dic["norst"])):
        sgas = list(dic["sgas"][j])
        rhog = list(dic["gas_den"][j])
        k_r = list(dic["gaskr"][j])
        rhow = list(dic["oil_den"][j])
        r_s = list(dic["rs"][j])
        krp = [k_rr > 0 for k_rr in k_r]
        krm = [k_rr <= 0 for k_rr in k_r]
        co2_g = [sga * rho * por for (sga, rho, por) in zip(sgas, rhog, dic["porva"])]
        co2_d = [
            rss * rho * (1.0 - sga) * por * GAS_DEN_REF / WAT_DEN_REF
            for (rss, rho, sga, por) in zip(r_s, rhow, sgas, dic["porva"])
        ]
        dic["ip1c"] = dic["pressure"][j][dic["sensor1"]] * 1e5  # Pa
        dic["ip2c"] = dic["pressure"][j][dic["sensor2"]] * 1e5
        dic["moba"] = sum((co2_g[i] * krp[i] for i in dic["boxa"]))
        dic["imma"] = sum((co2_g[i] * krm[i] for i in dic["boxa"]))
        dic["dissa"] = sum((co2_d[i] for i in dic["boxa"]))
        dic["seala"] = sum(
            ((co2_g[i] + co2_d[i]) * dic["facie1"][i] for i in dic["boxa"])
        )
        dic["mobb"] = sum((co2_g[i] * krp[i] for i in dic["boxb"]))
        dic["immb"] = sum((co2_g[i] * krm[i] for i in dic["boxb"]))
        dic["dissb"] = sum((co2_d[i] for i in dic["boxb"]))
        dic["sealb"] = sum(
            ((co2_g[i] + co2_d[i]) * dic["facie1"][i] for i in dic["boxb"])
        )
        dic["sealtot"] = sum(((co2_g[i] + co2_d[i]) for i in dic["facie1ind"]))
        if dic["case"] != "spe11c":
            dic["m_c"] = sum(
                abs((r_s[i_x] - r_s[i]) * dic["dz"][i])
                + abs((r_s[i_z] - r_s[i]) * dic["dx"][i])
                for (i, i_x, i_z) in zip(dic["boxc"], dic["boxc_x"], dic["boxc_z"])
            )
        else:
            dic["m_c"] = sum(
                abs((r_s[i_x] - r_s[i]) * dic["dy"][i] * dic["dz"][i])
                + abs((r_s[i_y] - r_s[i]) * dic["dx"][i] * dic["dz"][i])
                + abs((r_s[i_z] - r_s[i]) * dic["dx"][i] * dic["dy"][i])
                for (i, i_x, i_y, i_z) in zip(
                    dic["boxc"], dic["boxc_x"], dic["boxc_y"], dic["boxc_z"]
                )
            )
        text.append(
            f"{dic['d_t']*tim:.3e},{dic['ip1c']:.3e},{dic['ip2c']:.3e}"
            f",{dic['moba']:.3e},{dic['imma']:.3e},{dic['dissa']:.3e}"
            f",{dic['seala']:.3e},{dic['mobb']:.3e},{dic['immb']:.3e}"
            f",{dic['dissb']:.3e},{dic['sealb']:.3e},{dic['m_c']:.3e}"
            f",{dic['sealtot']:.3e}"
        )
        if dic["case"] != "spe11a":
            dic["boundtot"] = sum(((co2_g[i] + co2_d[i]) for i in dic["boundariesind"]))
            text[-1] += f",{dic['boundtot']:.3e}"

    with open(
        f"{dic['where']}/{dic['case']}_time_series.csv", "w", encoding="utf8"
    ) as file:
        file.write("\n".join(text))


def dense_data(dic):
    """Write the quantities with the benchmark format"""
    dic["d_t"] = (
        dic["smspec"].end_time - dic["unrst"].dates[dic["t0_rst"]]
    ).total_seconds() / (len(dic["unrst"].dates) - 1 - dic["t0_rst"])
    d_s = round(dic["spatial_t"] / dic["d_t"])
    dic["d_t"] = round(
        (dic["smspec"].end_time - dic["unrst"].dates[dic["t1_rst"]]).total_seconds()
        / dic["spatial_t"]
    )
    for i, j, k in zip(["x", "y", "z"], dic["dims"], dic["nxyz"]):
        temp = np.linspace(0, j, k + 1)
        dic[f"ref{i}cent"] = 0.5 * (temp[1:] + temp[:-1])
        dic[f"sim{i}cent"] = [0] * dic["nocellst"]
    ind, dic["cell_ind"] = 0, [0] * dic["nocellsr"]
    for cell in dic["egrid"].cells():
        (
            dic["simxcent"][cell.global_index],
            dic["simycent"][cell.global_index],
            dic["simzcent"][cell.global_index],
        ) = (
            cell.coordinate[0],
            cell.coordinate[1],
            dic["dims"][2] - cell.coordinate[2],
        )
    # This is very slow, it needs to be fixed
    for k in np.flip(dic["refzcent"]):
        for j in dic["refycent"]:
            for i in dic["refxcent"]:
                dic["cell_ind"][ind] = pd.Series(
                    (dic["simxcent"] - i) ** 2
                    + (dic["simycent"] - j) ** 2
                    + (dic["simzcent"] - k) ** 2
                ).argmin()
                ind += 1
    names = ["pressure", "sgas", "xco2", "xh20", "gden", "wden", "tco2"]
    if dic["case"] != "spe11a":
        names += ["temp"]
    for i in range(dic["d_t"] + 1):
        t_n = i * d_s + dic["t1_rst"]
        for name in names:
            dic[f"{name}_array"] = np.array([np.nan] * dic["nocellst"])
        co2_g = [
            sga * rho * por
            for (sga, rho, por) in zip(
                dic["sgas"][t_n], dic["gas_den"][t_n], dic["porva"]
            )
        ]
        co2_d = [
            rss * rho * (1.0 - sga) * por * GAS_DEN_REF / WAT_DEN_REF
            for (rss, rho, sga, por) in zip(
                dic["rs"][t_n], dic["oil_den"][t_n], dic["sgas"][t_n], dic["porva"]
            )
        ]
        h2o_l = [
            (1 - sga) * rho * por
            for (sga, rho, por) in zip(
                dic["sgas"][t_n], dic["oil_den"][t_n], dic["porva"]
            )
        ]
        dic["pressure_array"][dic["actind"]] = dic["pressure"][t_n] * 1e5  # Pa
        dic["sgas_array"][dic["actind"]] = dic["sgas"][t_n]
        dic["gden_array"][dic["actind"]] = [
            den if gas > 0 else 0
            for den, gas in zip(dic["gas_den"][t_n], dic["sgas"][t_n])
        ]
        dic["wden_array"][dic["actind"]] = dic["oil_den"][t_n]
        dic["xco2_array"][dic["actind"]] = [
            co2 / (co2 + h2o) for (co2, h2o) in zip(co2_d, h2o_l)
        ]
        dic["xh20_array"][dic["actind"]] = 0 * dic["xco2_array"][dic["actind"]]
        dic["tco2_array"][dic["actind"]] = [
            co2d + co2g for (co2d, co2g) in zip(co2_d, co2_g)
        ]
        if dic["case"] != "spe11a":
            dic["temp_array"][dic["actind"]] = dic["temp"][t_n]
            h2o_v = [
                rvv * rho * sga * por * WAT_DEN_REF / GAS_DEN_REF
                for (rvv, rho, sga, por) in zip(
                    dic["rv"][t_n], dic["gas_den"][t_n], dic["sgas"][t_n], dic["porva"]
                )
            ]
            dic["xh20_array"][dic["actind"]] = [
                h2o / (h2o + co2) if (h2o + co2) > 0 else 0.0
                for (h2o, co2) in zip(h2o_v, co2_g)
            ]
        for name in names:
            dic[f"{name}_refg"] = [np.nan] * dic["nocellsr"]
            dic[f"{name}_refg"] = dic[f"{name}_array"][dic["cell_ind"]]
        write_dense_data(dic, i)


def write_dense_data(dic, i):
    """Map the quantities to the cells"""
    if dic["case"] == "spe11a":
        name_t = f"{i*round(dic['spatial_t']/3600)}h"
        text = [
            "# x [m], z [m], pressure [Pa], gas saturation [-], "
            + "mass fraction of CO2 in liquid [-], mass fraction of H20 in vapor [-], "
            + "phase mass density gas [kg/m3], phase mass density water [kg/m3], "
            + "total mass CO2 [kg]"
        ]
    elif dic["case"] == "spe11b":
        name_t = f"{i*round(dic['spatial_t']/SECONDS_IN_YEAR)}y"
        text = [
            "# x [m], z [m], pressure [Pa], gas saturation [-], "
            + "mass fraction of CO2 in liquid [-], mass fraction of H20 in vapor [-], "
            + "phase mass density gas [kg/m3], phase mass density water [kg/m3], "
            + "total mass CO2 [kg], temperature [C]"
        ]
    else:
        name_t = f"{i*round(dic['spatial_t']/SECONDS_IN_YEAR)}y"
        text = [
            "# x [m], y [m], z [m], pressure [Pa], gas saturation [-], "
            + "mass fraction of CO2 in liquid [-], mass fraction of H20 in vapor [-], "
            + "phase mass density gas [kg/m3], phase mass density water [kg/m3], "
            + "total mass CO2 [kg], temperature [C]"
        ]
    idz = 0
    for zcord in dic["refzcent"]:
        idxy = 0
        for ycord in dic["refycent"]:
            for xcord in dic["refxcent"]:
                idc = (
                    dic["nxyz"][0] * dic["nxyz"][1] * (dic["nxyz"][2] - idz - 1) + idxy
                )
                if dic["case"] == "spe11a":
                    text.append(
                        f"{xcord:.3e}, {zcord:.3e}, {dic['pressure_refg'][idc] :.3e}, "
                        + f"{dic['sgas_refg'][idc] :.3e}, {dic['xco2_refg'][idc] :.3e}, "
                        + f"{dic['xh20_refg'][idc] :.3e}, {dic['gden_refg'][idc] :.3e}, "
                        + f"{dic['wden_refg'][idc] :.3e}, {dic['tco2_refg'][idc] :.3e}"
                    )
                elif dic["case"] == "spe11b":
                    text.append(
                        f"{xcord:.3e}, {zcord:.3e}, {dic['pressure_refg'][idc] :.3e}, "
                        + f"{dic['sgas_refg'][idc] :.3e}, {dic['xco2_refg'][idc] :.3e}, "
                        + f"{dic['xh20_refg'][idc] :.3e}, {dic['gden_refg'][idc] :.3e}, "
                        + f"{dic['wden_refg'][idc] :.3e}, {dic['tco2_refg'][idc] :.3e}, "
                        + f"{dic['temp_refg'][idc] :.3e}"
                    )
                else:
                    text.append(
                        f"{xcord:.3e}, {ycord:.3e}, {zcord:.3e}, "
                        + f"{dic['pressure_refg'][idc] :.3e}, {dic['sgas_refg'][idc] :.3e}, "
                        + f"{dic['xco2_refg'][idc] :.3e}, {dic['xh20_refg'][idc] :.3e}, "
                        + f"{dic['gden_refg'][idc] :.3e}, {dic['wden_refg'][idc] :.3e}, "
                        + f"{dic['tco2_refg'][idc] :.3e}, {dic['temp_refg'][idc] :.3e}"
                    )
                idxy += 1
        idz += 1
    with open(
        f"{dic['where']}/{dic['case']}_spatial_map_{name_t}.csv", "w", encoding="utf8"
    ) as file:
        file.write("\n".join(text))


if __name__ == "__main__":
    main()
