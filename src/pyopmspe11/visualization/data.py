# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT
# pylint: disable=C0302,R0912,R0914,R0801,R0915,E1102,C0325,R0902,R0913,R0917,R0911

"""Script to write the benchmark data"""

import argparse
import csv
import sys
from io import StringIO
from dataclasses import dataclass
from shapely.geometry import Polygon
from alive_progress import alive_bar
from rtree import index
import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import interp1d
from opm.io.ecl import EclFile as OpmFile
from opm.io.ecl import EGrid as OpmGrid
from opm.io.ecl import ERst as OpmRestart
from opm.io.ecl import ESmry as OpmSummary

GAS_DEN_REF = 1.86843
WAT_DEN_REF = 998.108
SECONDS_IN_YEAR = 31536000
SGAS_THR = 0.097


@dataclass(slots=True)
class BenchmarkConfig:
    """Paths and benchmark settings"""

    outfol: str
    case: str
    mode: str
    lower: bool
    deckfol: str
    flowfol: str
    where: str
    nxyz: NDArray
    dims: list
    denset: NDArray
    sparset: float
    nocellsrepgrid: int


@dataclass(slots=True)
class SimulationData:
    """Simulation data"""

    simres: str
    unrst: OpmRestart
    init: OpmFile
    egrid: OpmGrid
    smspec: OpmSummary
    times: list
    timesumary: list
    timeini: float
    noskiprst: int
    norst: int
    porv: NDArray
    porva: NDArray
    actind: list
    immiscible: bool
    isothermal: bool
    cornpoint: bool
    nocellst: int
    nocellsa: int
    dof: int
    nocellsxz: int
    simdim: list


def main(argv=None) -> None:
    """Postprocess simulation results into benchmark CSVs"""
    parser = argparse.ArgumentParser(description="Main script to process the data")
    parser.add_argument("-p", "--path", default="output", help="Output folder")
    parser.add_argument("-d", "--deck", default="spe11b", help="Simulated case")
    parser.add_argument("-r", "--resolution", default="10,1,5", help="x,y,z elements")
    parser.add_argument(
        "-t",
        "--time",
        default="24",
        help="Dense output time(s): spe11a [h], spe11b/c [y]",
    )
    parser.add_argument(
        "-w",
        "--write",
        default="0.1",
        help="Sparse/performance interval: spe11a [h], spe11b/c [y]",
    )
    parser.add_argument(
        "-g",
        "--generate",
        default="sparse",
        help="dense, sparse, performance, performance-spatial or combinations",
    )
    parser.add_argument(
        "-n", "--neighbourhood", default="", help="Region: 'lower' or all"
    )
    parser.add_argument("-f", "--subfolders", default=1, help="Create subfolders")
    cmdargs = vars(parser.parse_known_args(argv)[0])
    cfg = build_config_from_args(cmdargs)
    sim = read_simulations(cfg)

    if cfg.mode in (
        "performance",
        "all",
        "dense_performance",
        "performance_sparse",
        "dense_performance_sparse",
    ):
        performance(cfg, sim)
    if cfg.mode in (
        "all",
        "sparse",
        "dense_sparse",
        "dense_performance_sparse",
        "performance_sparse",
    ):
        sparse_data(cfg, sim)
    if cfg.mode in (
        "all",
        "performance-spatial",
        "dense",
        "dense_performance",
        "dense_sparse",
        "dense_performance-spatial",
        "dense_performance_sparse",
    ):
        if isinstance(cfg.denset, float):
            dt = cfg.denset
            cfg.denset = [i * dt for i in range(int(np.floor(sim.times[-1] / dt)) + 1)]
        dense_data(cfg, sim)
    print(f"The csv files have been written to {cfg.where}")


def build_config_from_args(cmdargs: dict) -> BenchmarkConfig:
    """Prepare the methods according to the spe11x case"""
    outfol = cmdargs["path"].strip()
    case = cmdargs["deck"].strip()
    mode = cmdargs["generate"].strip()
    lower = bool(cmdargs["neighbourhood"].strip())
    if int(cmdargs["subfolders"]) == 1:
        deckfol = f"{outfol}/deck"
        flowfol = f"{outfol}/flow"
        where = f"{outfol}/data"
    else:
        deckfol = flowfol = where = outfol
    nxyz = np.genfromtxt(StringIO(cmdargs["resolution"]), delimiter=",", dtype=int)
    if case == "spe11a":
        denset = (
            np.genfromtxt(StringIO(cmdargs["time"]), delimiter=",", dtype=float) * 3600
        )
        sparset = float(round(float(cmdargs["write"].strip()) * 3600))
        dims = [2.8, 1.0, 1.2]
        nxyz[1] = 1
    else:
        denset = (
            np.genfromtxt(StringIO(cmdargs["time"]), delimiter=",", dtype=float)
            * SECONDS_IN_YEAR
        )
        sparset = float(cmdargs["write"].strip()) * SECONDS_IN_YEAR
        dims = [8400.0, 1.0, 1200.0]
    if case == "spe11c":
        dims[1] = 5000.0
    nocellsrepgrid = nxyz[0] * nxyz[1] * nxyz[2]
    return BenchmarkConfig(
        outfol=outfol,
        case=case,
        mode=mode,
        lower=lower,
        deckfol=deckfol,
        flowfol=flowfol,
        where=where,
        nxyz=nxyz,
        dims=dims,
        denset=denset,
        sparset=sparset,
        nocellsrepgrid=nocellsrepgrid,
    )


def read_simulations(cfg: BenchmarkConfig) -> SimulationData:
    """Use opm Python package to read the results"""
    simres = f"{cfg.flowfol}/{cfg.outfol.split('/')[-1].upper()}"
    unrst = OpmRestart(f"{simres}.UNRST")
    immiscible = unrst.count("RSW", 0) == 0
    isothermal = unrst.count("TEMP", 0) == 0
    dof = 2 if isothermal else 3
    time = []
    times: list[float] = []
    timeini = 0
    noskiprst = 0
    for i in range(len(unrst.report_steps)):
        t = 86400 * unrst["DOUBHEAD", i][0]
        time.append(t)
        if not times:
            if (not immiscible and np.max(unrst["RSW", i]) > 0) or (
                immiscible and np.max(unrst["SGAS", i]) > 0
            ):
                timeini = 86400 * unrst["DOUBHEAD", i - 1][0]
                noskiprst = i - 1
                times = [0, t - timeini]
        else:
            times.append(t - timeini)
    if not times:
        times = time
    init = OpmFile(f"{simres}.INIT")
    egrid = OpmGrid(f"{simres}.EGRID")
    smspec = OpmSummary(f"{simres}.SMSPEC")
    norst = len(unrst.report_steps)
    porv = np.array(init["PORV"])
    actind = [i for i, p in enumerate(porv) if p > 0]
    porva = np.array([p for p in porv if p > 0])
    nocellst = len(porv)
    nocellsa = egrid.active_cells
    timesumary = [0.0] + list(86400.0 * smspec["TIME"] - timeini)
    dims = egrid.dimension
    simdim = [dims[0], dims[1], dims[2]]
    nocellsxz = dims[0] * dims[2]
    cornpoint = porv[-1] == 0
    return SimulationData(
        simres=simres,
        unrst=unrst,
        immiscible=immiscible,
        isothermal=isothermal,
        dof=dof,
        times=times,
        timeini=timeini,
        noskiprst=noskiprst,
        init=init,
        egrid=egrid,
        smspec=smspec,
        norst=norst,
        porv=porv,
        actind=actind,
        porva=porva,
        nocellst=nocellst,
        nocellsa=nocellsa,
        timesumary=timesumary,
        simdim=simdim,
        nocellsxz=nocellsxz,
        cornpoint=cornpoint,
    )


def performance(cfg: BenchmarkConfig, sim: SimulationData) -> None:
    """Generate benchmark performance data"""
    perf = build_performance_data(cfg, sim)
    write_performance_csv(cfg, perf)


def read_infostep(cfg: BenchmarkConfig, sim: SimulationData) -> tuple[list, NDArray]:
    """Read INFOSTEP file"""
    infosteps = []
    with open(
        f"{cfg.flowfol}/{cfg.outfol.split('/')[-1].upper()}.INFOSTEP",
        "r",
        encoding="utf8",
    ) as file:
        reader = csv.reader(file)
        tags = next(reader)[0].strip().split()
        for row in reader:
            values = row[0].strip().split()
            if float(values[0]) >= (sim.timeini - cfg.sparset) / 86400.0:
                infosteps.append([float(val) for val in values])
    return tags, np.array(infosteps)


def build_performance_data(cfg: BenchmarkConfig, sim: SimulationData) -> dict:
    """Build performance CSV data"""
    tags, infosteps = read_infostep(cfg, sim)
    infotimes = infosteps[:, tags.index("Time(day)")] * 86400.0 - sim.timeini
    times_data = np.linspace(0, sim.times[-1], round(sim.times[-1] / cfg.sparset) + 1)
    time_offset = max(0, sim.noskiprst - 1)
    map_info = np.array(
        [time_offset + int(np.floor(time_val / cfg.sparset)) for time_val in infotimes]
    )
    tmp = [0]
    for i in range(len(infotimes) - 1):
        if infotimes[i] != infotimes[i + 1]:
            tmp.append(tmp[-1] + 1)
        else:
            tmp.append(tmp[-1])
    detail_info = np.array(tmp)
    times_det = np.array(
        [np.max(infotimes[detail_info == i]) for i in range(np.max(detail_info) + 1)]
    )
    metrics = extract_solver_metrics(infosteps, tags)
    cpu_times, map_summary, summary_times = compute_cpu_times(cfg, sim, times_det)
    fgmip_values = sim.smspec["FGMIP"]
    if sim.timeini == 0:
        summary_times = np.insert(summary_times, 0, 0.0)
        fgmip_values = np.insert(fgmip_values, 0, 0.0)
    interp_fgmip = interp1d(summary_times, fgmip_values, fill_value="extrapolate")
    return {
        "series": build_time_series(
            sim, times_data, metrics, map_info, map_summary, interp_fgmip, cpu_times
        ),
        "detailed": build_detailed_series(
            sim, metrics, detail_info, infotimes, interp_fgmip, cpu_times
        ),
    }


def extract_solver_metrics(infosteps: NDArray, tags: list) -> dict:
    """Extract solver metrics"""
    return {
        "fsteps": np.array(infosteps[:, tags.index("Conv")] == 0, dtype=float),
        "nres": infosteps[:, tags.index("Lins")],
        "tlsolve": infosteps[:, tags.index("LSolve")],
        "linit": infosteps[:, tags.index("LinIt")],
        "nlinit": infosteps[:, tags.index("NewtIt")],
        "tsteps": 86400.0
        * infosteps[:, tags.index("TStep(day)")]
        * infosteps[:, tags.index("Conv")],
    }


def compute_cpu_times(
    cfg: BenchmarkConfig, sim: SimulationData, times_det: NDArray
) -> tuple[list, NDArray, NDArray]:
    """Compute CPU times"""
    cpu = sim.smspec["TCPU"]
    summary_times = 86400.0 * sim.smspec["TIME"] - sim.timeini
    map_summary = np.array(
        [
            max(0, sim.noskiprst - 1) + int(np.floor(time_val / cfg.sparset))
            for time_val in times_det
        ]
    )
    # For some spe11a cases (e.g., benchmark/spe11a/r3_cp_1cmish_capmax2500Pa.toml), it
    # seems a bug in OPM Flow is triggered where these arrays have different sizes. This is
    # a temporal fix, removing it (hopefully) later when that bug is fixed in OPM Flow.
    if cfg.case == "spe11a" and len(times_det) != len(cpu):
        interp_cpu = interp1d(summary_times, cpu, fill_value="extrapolate")
        tmp = np.append(times_det[1:], summary_times[-1])
        cpu = interp_cpu(tmp)
    if sim.timeini > 0:
        cpu = cpu[-len(map_summary) - 1 :]
        cpu = cpu[1:] - cpu[:-1]
    else:
        cpu = cpu[-len(map_summary) :]
        cpu[1:] -= cpu[:-1]
    # Extend later the temporal fix to spe11b/c cases if this issue is observed
    assert len(map_summary) == len(
        cpu
    ), "Please raise an issue with 'spe11b/c cpu/infostep fix'"
    return cpu, map_summary, np.array(summary_times)


def build_time_series(
    sim: SimulationData,
    times_data: NDArray,
    metrics: dict,
    map_info: NDArray,
    map_summary: NDArray,
    interp_fgmip: interp1d,
    cpu: list,
) -> list:
    """Build time series rows"""
    header = (
        "# t [s], tstep [s], fsteps [-], mass [kg], dof [-], "
        + "nliter [-], nres [-], liniter [-], runtime [s], tlinsol [s]"
    )
    rows = [header]
    if sim.noskiprst == 0:
        rows.append(
            f"0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, "
            f"{sim.dof * sim.nocellsa:.3e}, 0.000e+00, 0.000e+00, "
            f"0.000e+00, 0.000e+00, 0.000e+00"
        )
        times_data = np.delete(times_data, 0)

    fsteps = metrics["fsteps"]
    nres = metrics["nres"]
    tlsolve = metrics["tlsolve"]
    linit = metrics["linit"]
    nlinit = metrics["nlinit"]
    tsteps = metrics["tsteps"]
    dof_val = sim.dof * sim.nocellsa
    nrows = len(times_data)

    freq = np.zeros(nrows, dtype=int)
    run, itd = 0, 0
    for j in range(nrows):
        if np.sum(cpu[map_summary == j]) == 0:
            run += 1
            freq[j] = run
        else:
            run = 0
            freq[j] = 0

    weig = np.ones(nrows, dtype=float)
    quan = 1.0
    for j in range(nrows - 1, -1, -1):
        if freq[j] > 0 and quan == 1.0:
            quan = freq[j] + 1.0
        elif freq[j] == 0:
            weig[j] = quan
            quan = 1.0
            continue
        weig[j] = quan

    max_block = map_info.max() + 1
    sum_tstep = np.zeros(max_block)
    cnt_tstep = np.zeros(max_block, dtype=int)
    sum_fsteps = np.zeros(max_block)
    sum_nlinit = np.zeros(max_block)
    sum_nres = np.zeros(max_block)
    sum_linit = np.zeros(max_block)
    sum_tlinsol = np.zeros(max_block)

    for b in range(max_block):
        ind = map_info == b
        if np.any(ind):
            sum_tstep[b] = np.sum(tsteps[ind])
            cnt_tstep[b] = np.sum(ind)
            sum_fsteps[b] = np.sum(fsteps[ind])
            sum_nlinit[b] = np.sum(nlinit[ind])
            sum_nres[b] = np.sum(nres[ind])
            sum_linit[b] = np.sum(linit[ind])
            sum_tlinsol[b] = np.sum(tlsolve[ind])

    cur_block = None
    cur_tstep = 0.0

    for j, time_val in enumerate(times_data):
        if freq[j] == 0:
            cur_block = j
            itd = map_summary == j
            if cnt_tstep[j] > 0:
                cur_tstep = sum_tstep[j] / cnt_tstep[j]
            else:
                cur_tstep = 0.0
        w = weig[j]
        rows.append(
            f"{time_val:.3e}, "
            f"{cur_tstep / w:.3e}, "
            f"{sum_fsteps[cur_block] / w:.3e}, "
            f"{interp_fgmip(time_val):.3e}, "
            f"{dof_val:.3e}, "
            f"{sum_nlinit[cur_block] / w:.3e}, "
            f"{sum_nres[cur_block] / w:.3e}, "
            f"{sum_linit[cur_block] / w:.3e}, "
            f"{np.sum(cpu[itd]) / w:.3e}, "
            f"{sum_tlinsol[cur_block] / w:.3e}"
        )

    return rows


def build_detailed_series(
    sim: SimulationData,
    metrics: dict,
    detail_info: NDArray,
    infotimes: NDArray,
    interp_fgmip: interp1d,
    cpu: list,
) -> list:
    """Build detailed time series rows"""
    header = (
        "# t [s], tstep [s], fsteps [-], mass [kg], dof [-], nliter [-], "
        + "nres [-], liniter [-], runtime [s], tlinsol [s]"
    )
    rows = [header]
    for detail_index in range(np.max(detail_info) + 1):
        mask = detail_info == detail_index
        time_val = np.max(infotimes[mask])
        if time_val >= 0:
            rows.append(
                f"{time_val:.3e}, {np.max(metrics['tsteps'][mask]):.3e}, "
                f"{np.sum(metrics['fsteps'][mask]):.3e}, "
                f"{interp_fgmip(time_val):.3e}, {sim.dof*sim.nocellsa:.3e}, "
                f"{np.sum(metrics['nlinit'][mask]):.3e}, {np.sum(metrics['nres'][mask]):.3e}, "
                f"{np.sum(metrics['linit'][mask]):.3e}, {cpu[detail_index]:.3e}, "
                f"{np.sum(metrics['tlsolve'][mask]):.3e}"
            )
    return rows


def write_performance_csv(cfg: BenchmarkConfig, perf: dict) -> None:
    """Write performance CSV files"""
    with open(
        f"{cfg.where}/{cfg.case}_performance_time_series.csv", "w", encoding="utf8"
    ) as file:
        file.write("\n".join(perf["series"]))
    with open(
        f"{cfg.where}/{cfg.case}_performance_time_series_detailed.csv",
        "w",
        encoding="utf8",
    ) as file:
        file.write("\n".join(perf["detailed"]))


def sparse_data(cfg: BenchmarkConfig, sim: SimulationData) -> None:
    """Generate sparse benchmark data"""
    sparse = build_sparse_data(cfg, sim)
    write_sparse_csv(cfg, sparse)


def build_sparse_data(cfg: BenchmarkConfig, sim: SimulationData) -> dict:
    """Build sparse benchmark data"""
    times_data = np.linspace(0, sim.times[-1], round(sim.times[-1] / cfg.sparset) + 1)
    fipnum = list(sim.init["FIPNUM"])
    dx = np.array(sim.init["DX"])
    dy = np.array(sim.init["DY"])
    dz = np.array(sim.init["DZ"])
    fip_groups = get_fip_groups(cfg)
    summary_data = build_summary_data(cfg, sim, fipnum, fip_groups)
    m_c = (
        compute_m_c(cfg, sim, fipnum, dx, dy, dz)
        if not sim.immiscible
        else [0.0] * (sim.norst - sim.noskiprst - 1)
    )
    interpolated = interpolate_sparse(times_data, sim, summary_data, m_c)
    return interpolated


def get_fip_groups(cfg: BenchmarkConfig) -> dict:
    """Define FIP groups"""
    if cfg.lower:
        result = {
            "diss_a": [2, 4, 8],
            "seal_a": [],
            "diss_b": [],
            "seal_b": [],
            "bound": [],
        }
    else:
        result = {
            "diss_a": [2, 4, 5, 8],
            "seal_a": [5, 8],
            "diss_b": [3, 6],
            "seal_b": [6],
            "bound": [],
        }
    if cfg.case != "spe11a":
        result["bound"] = [11]
    if cfg.case == "spe11c":
        result["diss_a"] += [13, 17]
        result["bound"] += [13, 17]
        if not cfg.lower:
            result["diss_a"] += [14]
            result["seal_a"] += [14]
            result["diss_b"] += [15, 16]
            result["seal_b"] += [16]
            result["bound"] += [14, 15, 16]
    return result


def build_summary_data(
    cfg: BenchmarkConfig, sim: SimulationData, fipnum: list, groups: dict
) -> dict:
    """Build sparse summary time series"""
    zero_series = 0.0 * sim.smspec["TIME"]
    pop1, pop2 = extract_boundary_pressures(cfg, sim, fipnum)
    result = {
        "pop1": pop1,
        "pop2": pop2,
        "moba": zero_series.copy(),
        "imma": zero_series.copy(),
        "dissa": zero_series.copy(),
        "seala": zero_series.copy(),
        "mobb": zero_series.copy(),
        "immb": zero_series.copy(),
        "dissb": zero_series.copy(),
        "sealb": zero_series.copy(),
        "sealt": zero_series.copy(),
        "boundtot": zero_series.copy(),
    }
    for i in groups["diss_a"]:
        result["moba"] += sim.smspec[f"RGKMO:{i}"]
        result["imma"] += sim.smspec[f"RGKTR:{i}"]
        result["dissa"] += sim.smspec[f"RGMDS:{i}"]
    for i in groups["seal_a"]:
        result["seala"] += sim.smspec[f"RGMDS:{i}"]
        result["seala"] += sim.smspec[f"RGKMO:{i}"]
        result["seala"] += sim.smspec[f"RGKTR:{i}"]
    for i in groups["diss_b"]:
        result["mobb"] += sim.smspec[f"RGKMO:{i}"]
        result["immb"] += sim.smspec[f"RGKTR:{i}"]
        result["dissb"] += sim.smspec[f"RGMDS:{i}"]
    for i in groups["seal_b"]:
        result["sealb"] += sim.smspec[f"RGMDS:{i}"]
        result["sealb"] += sim.smspec[f"RGKMO:{i}"]
        result["sealb"] += sim.smspec[f"RGKTR:{i}"]
    result["sealt"] = result["seala"] + result["sealb"]
    for name in ("RGMDS", "RGKMO", "RGKTR"):
        if not cfg.lower:
            result["sealt"] += sim.smspec[f"{name}:7"]
            result["sealt"] += sim.smspec[f"{name}:9"]
        else:
            if 7 in fipnum:
                result["sealt"] += sim.smspec[f"{name}:7"]
    for i in groups["bound"]:
        if i in fipnum:
            result["boundtot"] += sim.smspec[f"RGMDS:{i}"]
            result["boundtot"] += sim.smspec[f"RGKMO:{i}"]
            result["boundtot"] += sim.smspec[f"RGKTR:{i}"]
    return result


def extract_boundary_pressures(
    cfg: BenchmarkConfig, sim: SimulationData, fipnum: list
) -> tuple[list, list]:
    """Extract boundary pressures"""
    pressure = sim.unrst["PRESSURE", 0]
    pcgw = sim.unrst["PCGW", 0]
    index_pop1 = fipnum.index(8)
    pop1_value = (pressure[index_pop1] - pcgw[index_pop1]) * 1.0e5
    summary_keys = [key for key in sim.smspec.keys() if key.startswith("BWPR")]
    summary_keys.sort()
    pop1 = [pop1_value] + list(sim.smspec[summary_keys[0]] * 1.0e5)
    if cfg.lower:
        pop2 = pop1
    else:
        index_pop2 = fipnum.index(9)
        pop2_value = (pressure[index_pop2] - pcgw[index_pop2]) * 1.0e5
        pop2 = [pop2_value] + list(sim.smspec[summary_keys[1]] * 1.0e5)
    return pop1, pop2


def compute_m_c(
    cfg: BenchmarkConfig,
    sim: SimulationData,
    fipnum: list,
    dx: NDArray,
    dy: NDArray,
    dz: NDArray,
) -> list:
    """Compute Box C variation"""
    box_mask = np.isin(fipnum, (4, 12, 17, 18))
    box_x = np.roll(box_mask, 1)
    box_y = np.roll(box_mask, -sim.simdim[0])
    box_z = np.roll(box_mask, -sim.simdim[0] * sim.simdim[1])
    den_ratio = WAT_DEN_REF / GAS_DEN_REF
    dx_b = dx[box_mask]
    dy_b = dy[box_mask]
    dz_b = dz[box_mask]
    values = []
    for step in range(sim.noskiprst + 1, sim.norst):
        rss = np.array(sim.unrst["RSW", step])
        max_sat = np.array(sim.unrst["RSWSAT", step])
        xcw = rss / (rss + den_ratio)
        xcw /= max_sat / (max_sat + den_ratio)
        dxv = np.abs(xcw[box_x] - xcw[box_mask])
        dzv = np.abs(xcw[box_z] - xcw[box_mask])
        if cfg.case != "spe11c":
            values.append(np.sum(dxv * dz_b + dzv * dx_b))
        else:
            dyv = np.abs(xcw[box_y] - xcw[box_mask])
            values.append(
                np.sum(dxv * dy_b * dz_b + dyv * dx_b * dz_b + dzv * dx_b * dy_b)
            )
    return values


def interpolate_sparse(
    times_data: NDArray, sim: SimulationData, summary: dict, m_c: list
) -> dict:
    """Interpolate sparse outputs"""
    result = {}
    tsim = sim.timesumary
    tlen = len(tsim)
    for key, values in summary.items():
        if isinstance(values, float):
            series = [0.0] * (tlen - 1)
        else:
            series = list(values)
        if len(series) == tlen - 1:
            series = [0.0] + series
        interp = interp1d(tsim, series, fill_value="extrapolate")
        result[key] = interp(times_data)
    interp_mc = interp1d(sim.times, [0.0] + list(m_c), fill_value="extrapolate")
    result["m_c"] = interp_mc(times_data)
    result["times"] = times_data
    return result


def write_sparse_csv(cfg: BenchmarkConfig, sparse: dict) -> None:
    """Write sparse CSV"""
    header = (
        "# t [s], p1 [Pa], p2 [Pa], mobA [kg], immA [kg], dissA [kg], sealA [kg], "
        + "mobB [kg], immB [kg], dissB [kg], sealB [kg], MC [m], sealTot [kg]"
    )
    times = sparse["times"]
    pop1 = sparse["pop1"]
    pop2 = sparse["pop2"]
    moba = sparse["moba"]
    imma = sparse["imma"]
    dissa = sparse["dissa"]
    seala = sparse["seala"]
    mobb = sparse["mobb"]
    immb = sparse["immb"]
    dissb = sparse["dissb"]
    sealb = sparse["sealb"]
    mc = sparse["m_c"]
    sealt = sparse["sealt"]
    path = f"{cfg.where}/{cfg.case}_time_series.csv"
    with open(path, "w", encoding="utf8") as file:
        if cfg.case == "spe11a":
            file.write(header)
            for i, time_val in enumerate(times):
                file.write(
                    f"\n{time_val:.3e}, {pop1[i]:.5e}, {pop2[i]:.5e}, "
                    f"{moba[i]:.3e}, {imma[i]:.3e}, {dissa[i]:.3e}, "
                    f"{seala[i]:.3e}, {mobb[i]:.3e}, {immb[i]:.3e}, "
                    f"{dissb[i]:.3e}, {sealb[i]:.3e}, {mc[i]:.3e}, "
                    f"{sealt[i]:.3e}"
                )
        else:
            file.write(header + ", boundTot [kg]")
            boundtot = sparse["boundtot"]
            for i, time_val in enumerate(times):
                file.write(
                    f"\n{time_val:.4e}, {pop1[i]:.3e}, {pop2[i]:.3e}, "
                    f"{moba[i]:.3e}, {imma[i]:.3e}, {dissa[i]:.3e}, "
                    f"{seala[i]:.3e}, {mobb[i]:.3e}, {immb[i]:.3e}, "
                    f"{dissb[i]:.3e}, {sealb[i]:.3e}, {mc[i]:.3e}, "
                    f"{sealt[i]:.3e}, {boundtot[i]:.3e}"
                )


def dense_data(cfg: BenchmarkConfig, sim: SimulationData) -> None:
    """Generate dense benchmark data"""
    rstno, nrstno, refgrid, mapping, actindr = build_dense_static(cfg, sim)
    if cfg.mode == "all" or cfg.mode[:5] == "dense":
        for step_index, restart in enumerate(rstno):
            if step_index + 1 < nrstno:
                print(f"Processing dense data {step_index+1} out of {nrstno}", end="\r")
            else:
                print(f"Processing dense data {step_index+1} out of {nrstno}")
            restart_index = restart + sim.noskiprst
            dense_step = build_dense_step(cfg, sim, mapping, restart_index, actindr)
            write_dense_csv(cfg, sim, refgrid, dense_step, step_index)
    if cfg.mode in ("all", "performance-spatial", "dense_performance-spatial"):
        handle_performance_spatial(cfg, sim, rstno, refgrid, mapping, actindr)


def can_use_fast_dense_mapping(
    cfg: BenchmarkConfig, sim: SimulationData, dx: NDArray, dz: NDArray
) -> bool:
    """Check if fast dense mapping can be used"""
    if cfg.lower:
        return False
    if np.min(dz) != np.max(dz):
        return False
    if cfg.nxyz[2] == sim.simdim[2] and cfg.nxyz[0] == sim.simdim[0]:
        return np.min(dx) == np.max(dx)
    if cfg.nxyz[2] == sim.simdim[2] and cfg.nxyz[0] == sim.simdim[0] - 2:
        return np.min(dx[2:-2]) == np.max(dx[2:-2])
    if sim.simdim[2] % cfg.nxyz[2] == 0 and sim.simdim[0] % cfg.nxyz[0] == 0:
        return np.min(dx) == np.max(dx)
    if sim.simdim[2] % cfg.nxyz[2] == 0 and (sim.simdim[0] - 2) % cfg.nxyz[0] == 0:
        return np.min(dx[2:-2]) == np.max(dx[2:-2])
    return False


def build_general_dense_mapping(
    cfg: BenchmarkConfig,
    sim: SimulationData,
    refgrid: tuple[NDArray, NDArray, NDArray, NDArray, NDArray, NDArray],
    geometry: tuple[NDArray, NDArray, list],
) -> tuple[list[list[list[int | float]]], NDArray]:
    """General sim-to-report mapping using geometry"""
    refxvert, _, refzvert, refxcent, _, refzcent = refgrid
    simxcent, simzcent, simpoly = geometry
    cell_ind: list[list[list[int | float]]] = [[] for _ in range(sim.nocellsxz)]
    cell_cent = np.zeros(cfg.nocellsrepgrid, dtype=float)
    nrefx = cfg.nxyz[0]
    nrefz = cfg.nxyz[2]
    refpoly = [None] * (nrefx * nrefz)
    refxgrid = np.zeros(nrefx * nrefz)
    refzgrid = np.zeros(nrefx * nrefz)
    idx = index.Index()
    rid = 0
    for kz, zcen in enumerate(refzcent):
        zv0 = refzvert[kz]
        zv1 = refzvert[kz + 1]
        for ix, xcen in enumerate(refxcent):
            xv0 = refxvert[ix]
            xv1 = refxvert[ix + 1]
            refxgrid[rid] = xcen
            refzgrid[rid] = zcen
            poly = Polygon(((xv0, zv0), (xv1, zv0), (xv1, zv1), (xv0, zv1)))
            refpoly[rid] = poly
            idx.insert(rid, poly.bounds)
            rid += 1
    print("Processing polygon intersections between simulation and reporting grids")
    with alive_bar(len(simpoly), bar="fish") as bar_animation:
        for sim_cell, poly_s in enumerate(simpoly):
            bar_animation()
            if poly_s.area > 0.0:
                area_s = poly_s.area
                for tgt in idx.intersection(poly_s.bounds):
                    a = poly_s.intersection(refpoly[tgt]).area
                    if a > 0.0:
                        cell_ind[sim_cell].append([tgt, a / area_s])
            else:
                cell_ind[sim_cell] = cell_ind[sim_cell - 1]
    sx = np.asarray(simxcent)
    sz = np.asarray(simzcent)
    print("Finding the cell indices between simulation and reporting grids")
    with alive_bar(len(refxgrid), bar="fish") as bar_animation:
        for rep, (xc, zc) in enumerate(zip(refxgrid, refzgrid)):
            bar_animation()
            cell_cent[rep] = np.argmin(np.abs(sx - xc) + np.abs(sz - zc))
    return cell_ind, cell_cent


def build_dense_static(cfg: BenchmarkConfig, sim: SimulationData) -> tuple[
    list,
    int,
    tuple[NDArray, NDArray, NDArray, NDArray, NDArray, NDArray],
    tuple[list[list[list[int | float]]], NDArray],
    NDArray,
]:
    """Build static dense grid data"""
    rstno, nrstno = build_dense_schedule(cfg, sim)
    refgrid = build_dense_reference_grid(cfg)
    simxcent, simycent, simzcent, simpoly, satnum = extract_sim_geometry(cfg, sim)
    dx = np.array(sim.init["DX"])
    dz = np.array(sim.init["DZ"])
    if can_use_fast_dense_mapping(cfg, sim, dx, dz):
        cell_ind, cell_cent = build_fast_dense_mapping(cfg, sim, dx, dz)
    else:
        cell_ind, cell_cent = build_general_dense_mapping(
            cfg, sim, refgrid, (simxcent, simzcent, simpoly)
        )
    cell_ind, cell_cent, actindr = handle_post_mapping(
        cfg, sim, refgrid, (cell_ind, cell_cent), simycent, satnum
    )
    return rstno, nrstno, refgrid, (cell_ind, cell_cent), actindr


def build_dense_schedule(cfg: BenchmarkConfig, sim: SimulationData) -> tuple[list, int]:
    """Compute restart indices for dense output"""
    rstno = [sim.times.index(time_val) for time_val in cfg.denset]
    return rstno, len(rstno)


def build_dense_reference_grid(
    cfg: BenchmarkConfig,
) -> tuple[NDArray, NDArray, NDArray, NDArray, NDArray, NDArray]:
    """Build reporting grid coordinates"""
    refxvert = np.linspace(0, cfg.dims[0], cfg.nxyz[0] + 1)
    refyvert = np.linspace(0, cfg.dims[1], cfg.nxyz[1] + 1)
    refzvert = np.linspace(0, cfg.dims[2], cfg.nxyz[2] + 1)
    refxcent = 0.5 * (refxvert[1:] + refxvert[:-1])
    refycent = 0.5 * (refyvert[1:] + refyvert[:-1])
    refzcent = 0.5 * (refzvert[1:] + refzvert[:-1])
    return refxvert, refyvert, refzvert, refxcent, refycent, refzcent


def extract_sim_geometry(
    cfg: BenchmarkConfig, sim: SimulationData
) -> tuple[NDArray, NDArray, NDArray, list, NDArray]:
    """Extract simulation geometry"""
    simxcent = np.zeros(sim.nocellsxz)
    simzcent = np.zeros(sim.nocellsxz)
    simycent = np.zeros(sim.simdim[1])
    simpoly = [None] * sim.nocellsxz
    satnum = np.array(sim.init["SATNUM"])
    z_0 = 155.04166666666666 if cfg.case == "spe11c" else 0.0
    nx, ny, nz = sim.simdim
    dims_z = cfg.dims[2]
    for j in range(nz):
        for i in range(nx):
            n = i + (nz - j - 1) * nx
            xyz = sim.egrid.xyz_from_ijk(i, 0, nz - j - 1)
            poly = Polygon(
                [
                    [xyz[0][0], dims_z - (xyz[2][0] - z_0)],
                    [xyz[0][1], dims_z - (xyz[2][1] - z_0)],
                    [xyz[0][5], dims_z - (xyz[2][5] - z_0)],
                    [xyz[0][4], dims_z - (xyz[2][4] - z_0)],
                ]
            )
            simpoly[n] = poly
            xcen, zcen = (float(v) for v in poly.centroid.wkt[7:-1].split())
            simxcent[n], simzcent[n] = xcen, zcen
    for j in range(ny):
        xyz = sim.egrid.xyz_from_ijk(0, j, 0)
        simycent[j] = 0.5 * (xyz[1][2] - xyz[1][1]) + xyz[1][1]
    if cfg.lower and sim.cornpoint:
        nx = sim.simdim[0]

        simzcent = np.array(simzcent)
        simxcent = np.insert(simxcent, 0, simxcent[:nx])
        simzcent = np.insert(simzcent, 0, simzcent[:nx] + 1e-4)
    return simxcent, simycent, simzcent, simpoly, satnum


def build_fast_dense_mapping(
    cfg: BenchmarkConfig, sim: SimulationData, dx: NDArray, dz: NDArray
) -> tuple[list[list[list[int | float]]], NDArray]:
    """Fast sim-to-report mapping for aligned grids"""
    cell_ind: list[list[list[int | float]]] = [[] for _ in range(sim.nocellsxz)]
    cell_indc = np.zeros(sim.nocellsxz, dtype=int)
    cell_cent = np.zeros(cfg.nocellsrepgrid, dtype=float)
    iszunif = (np.min(dz) == np.max(dz)) and not cfg.lower

    if (
        iszunif
        and cfg.nxyz[2] == sim.simdim[2]
        and np.min(dx) == np.max(dx)
        and cfg.nxyz[0] == sim.simdim[0]
    ):
        for layer in range(cfg.nxyz[2]):
            cell_cent[
                (cfg.nxyz[2] - layer - 1)
                * cfg.nxyz[0] : (cfg.nxyz[2] - layer)
                * cfg.nxyz[0]
            ] = range(layer * cfg.nxyz[0], (layer + 1) * cfg.nxyz[0])
        cell_indc[:] = cell_cent
        for i, value in enumerate(cell_indc):
            cell_ind[i] = [[int(value), 1.0]]

    elif (
        iszunif
        and cfg.nxyz[2] == sim.simdim[2]
        and np.min(dx[2:-2]) == np.max(dx[2:-2])
        and cfg.nxyz[0] == sim.simdim[0] - 2
    ):
        for layer in range(cfg.nxyz[2]):
            base = (sim.simdim[2] - layer - 1) * sim.simdim[0]
            cell_indc[base] = layer * cfg.nxyz[0]
            cell_indc[base + 1 : base + sim.simdim[0] - 1] = range(
                layer * cfg.nxyz[0], (layer + 1) * cfg.nxyz[0] + 1
            )
            cell_indc[base + sim.simdim[0] - 1] = (layer + 1) * cfg.nxyz[0] - 1
            cell_cent[
                (cfg.nxyz[2] - layer - 1)
                * cfg.nxyz[0] : (cfg.nxyz[2] - layer)
                * cfg.nxyz[0]
            ] = [
                value + 2 * layer
                for value in range(
                    layer * cfg.nxyz[0] + 1, (layer + 1) * cfg.nxyz[0] + 1
                )
            ]
        for i, value in enumerate(cell_indc):
            cell_ind[i] = [[int(value), 1.0]]

    elif (
        iszunif
        and sim.simdim[2] % cfg.nxyz[2] == 0
        and np.min(dx) == np.max(dx)
        and sim.simdim[0] % cfg.nxyz[0] == 0
    ):
        x_repeat = sim.simdim[0] // cfg.nxyz[0]
        z_repeat = sim.simdim[2] // cfg.nxyz[2]
        for layer in range(cfg.nxyz[2]):
            cell_cent[
                (cfg.nxyz[2] - layer - 1)
                * cfg.nxyz[0] : (cfg.nxyz[2] - layer)
                * cfg.nxyz[0]
            ] = [
                value * x_repeat
                + (x_repeat / 2 - 1)
                + (z_repeat / 2 - 1) * sim.simdim[0]
                + layer * (z_repeat - 1) * sim.simdim[0]
                for value in range(layer * cfg.nxyz[0], (layer + 1) * cfg.nxyz[0])
            ]
        for layer in range(cfg.nxyz[2]):
            for zloc in range(z_repeat):
                for xloc in range(cfg.nxyz[0]):
                    start = (
                        sim.simdim[2] - (layer * z_repeat + zloc) - 1
                    ) * sim.simdim[0] + xloc * x_repeat
                    cell_indc[start : start + x_repeat] = [
                        xloc + layer * cfg.nxyz[0]
                    ] * x_repeat
        for i, value in enumerate(cell_indc):
            cell_ind[i] = [[int(value), 1.0]]

    elif (
        iszunif
        and sim.simdim[2] % cfg.nxyz[2] == 0
        and np.min(dx[2:-2]) == np.max(dx[2:-2])
        and (sim.simdim[0] - 2) % cfg.nxyz[0] == 0
    ):
        x_repeat = (sim.simdim[0] - 2) // cfg.nxyz[0]
        z_repeat = sim.simdim[2] // cfg.nxyz[2]
        for layer in range(cfg.nxyz[2]):
            cell_cent[
                (cfg.nxyz[2] - layer - 1)
                * cfg.nxyz[0] : (cfg.nxyz[2] - layer)
                * cfg.nxyz[0]
            ] = [
                value * x_repeat
                + (x_repeat / 2 - 1)
                + (z_repeat / 2 - 1) * sim.simdim[0]
                + layer * (z_repeat - 1) * sim.simdim[0]
                + 2 * layer
                + 1
                for value in range(layer * cfg.nxyz[0], (layer + 1) * cfg.nxyz[0])
            ]
        for layer in range(cfg.nxyz[2]):
            for zloc in range(z_repeat):
                for xloc in range(cfg.nxyz[0]):
                    base = (sim.simdim[2] - (layer * z_repeat + zloc) - 1) * sim.simdim[
                        0
                    ] + xloc * x_repeat
                    cell_indc[base] = layer * cfg.nxyz[0] + xloc
                    cell_indc[base + 1 : base + x_repeat + 1] = [
                        layer * cfg.nxyz[0] + xloc
                    ] * x_repeat
                    cell_indc[base + x_repeat + 1] = layer * cfg.nxyz[0] + xloc
        for i, value in enumerate(cell_indc):
            cell_ind[i] = [[int(value), 1.0]]

    return cell_ind, cell_cent


def handle_post_mapping(
    cfg: BenchmarkConfig,
    sim: SimulationData,
    refgrid: tuple[NDArray, NDArray, NDArray, NDArray, NDArray, NDArray],
    mapping: tuple[list[list[list[int | float]]], NDArray],
    simycent: NDArray,
    satnum: NDArray,
) -> tuple[list[list[list[int | float]]], NDArray, NDArray]:
    """Post-process mapping"""
    cell_ind, cell_cent = mapping
    actindr = np.empty(0)
    if np.max(satnum) < 7 and cfg.case == "spe11a":
        actindr = handle_inactive_mapping(cfg, sim, cell_ind)
    if cfg.case == "spe11c":
        cell_cent = handle_yaxis_mapping_intensive(
            cfg, sim, cell_cent, refgrid[4], simycent
        )
        cell_ind = handle_yaxis_mapping_extensive(
            cfg, sim, cell_ind, simycent, refgrid[1]
        )
    return cell_ind, cell_cent, actindr


def build_dense_step(
    cfg: BenchmarkConfig,
    sim: SimulationData,
    mapping: tuple[list[list[list[int | float]]], NDArray],
    restart_index: int,
    actindr: NDArray,
) -> dict:
    """Build dense data for one restart"""
    names = ["pressure", "sgas", "xco2", "xh20", "gden", "wden", "tco2"]
    if not sim.isothermal:
        names = ["temp"] + names
    arrays = generate_arrays(cfg, sim, names, restart_index, actindr)
    map_to_report_grid(sim, mapping, arrays)
    return arrays


def write_dense_csv(
    cfg: BenchmarkConfig,
    sim: SimulationData,
    refgrid: tuple[NDArray, NDArray, NDArray, NDArray, NDArray, NDArray],
    dense_step: dict,
    step_index: int,
) -> None:
    """Write dense spatial CSV"""
    name_t, text = get_header(cfg, sim, step_index)
    if cfg.lower:
        dense_step["tco2_refg"][np.isnan(dense_step["sgas_refg"])] = np.nan
    _, _, _, refxcent, refycent, refzcent = refgrid
    nx, ny, nz = cfg.nxyz
    p_arr = dense_step["pressure_refg"]
    s_arr = dense_step["sgas_refg"]
    xco2_arr = dense_step["xco2_refg"]
    xh20_arr = dense_step["xh20_refg"]
    gden_arr = dense_step["gden_refg"]
    wden_arr = dense_step["wden_refg"]
    tco2_arr = dense_step["tco2_refg"]
    temp_arr = dense_step.get("temp_refg")
    path = f"{cfg.where}/{cfg.case}_spatial_map_{name_t}.csv"
    with open(path, "w", encoding="utf8") as file:
        file.write("\n".join(text))
        idz = 0
        for zcord in refzcent:
            idxy = 0
            for ycord in refycent:
                for xcord in refxcent:
                    cell = -nx * ny * (nz - idz) + idxy
                    p = p_arr[cell]
                    if np.isnan(p):
                        co2v = tco2_arr[cell]
                        co2 = "n/a" if np.isnan(co2v) else f"{co2v:.3e}"
                        if cfg.case == "spe11a":
                            row = (
                                f"{xcord:.3e}, {zcord:.3e}, n/a, n/a, n/a, n/a, n/a, n/a, {co2}"
                                + (", n/a" if not sim.isothermal else "")
                            )
                        elif cfg.case == "spe11b":
                            row = (
                                f"{xcord:.3e}, {zcord:.3e}, n/a, n/a, n/a, n/a, n/a, n/a, "
                                f"{co2}, n/a"
                            )
                        else:
                            row = (
                                f"{xcord:.3e}, {ycord:.3e}, {zcord:.3e}, n/a, n/a, n/a, n/a, "
                                f"n/a, n/a, {co2}, n/a"
                            )
                    else:
                        pf = f"{p:.3e}"
                        sf = f"{s_arr[cell]:.3e}"
                        gf = f"{gden_arr[cell]:.3e}"
                        wf = f"{wden_arr[cell]:.3e}"
                        co2 = f"{tco2_arr[cell]:.3e}"
                        if sim.immiscible:
                            xf = hf = "n/a"
                        else:
                            xf = f"{xco2_arr[cell]:.3e}"
                            hf = f"{xh20_arr[cell]:.3e}"
                        tf = f"{temp_arr[cell]:.3e}" if temp_arr is not None else None
                        if cfg.case == "spe11a":
                            row = (
                                f"{xcord:.3e}, {zcord:.3e}, {pf}, {sf}, {xf}, {hf}, "
                                f"{gf}, {wf}, {co2}" + (f", {tf}" if tf else "")
                            )
                        elif cfg.case == "spe11b":
                            row = (
                                f"{xcord:.3e}, {zcord:.3e}, {pf}, {sf}, {xf}, {hf}, {gf}, "
                                f"{wf}, {co2}, {tf}"
                            )
                        else:
                            row = (
                                f"{xcord:.3e}, {ycord:.3e}, {zcord:.3e}, {pf}, {sf}, {xf}, "
                                f"{hf}, {gf}, {wf}, {co2}, {tf}"
                            )
                    file.write("\n" + row)
                    idxy += 1
            idz += 1


def handle_yaxis_mapping_extensive(
    cfg: BenchmarkConfig,
    sim: SimulationData,
    cell_ind: list[list[list[int | float]]],
    simycent: NDArray,
    refyvert: NDArray,
) -> list[list[list[int | float]]]:
    """Extend indices for y direction (extensive)"""
    simyvert = [0.0]
    for yval in simycent:
        simyvert.append(simyvert[-1] + 2.0 * (yval - simyvert[-1]))
    weights = []
    indy = []
    ind = 0
    for y_i, y_f in zip(simyvert[:-1], simyvert[1:]):
        if refyvert[ind + 1] <= y_i:
            ind += 1
        if refyvert[ind] <= y_i and y_f <= refyvert[ind + 1]:
            indy.append(ind)
            weights.append([1.0])
        else:
            w0 = (refyvert[ind + 1] - y_i) / (y_f - y_i)
            w1 = (y_f - refyvert[ind + 1]) / (y_f - y_i)
            indy.append(ind)
            weights.append([w0, w1])
            ind += 1
    nx, ny, nz = sim.simdim
    nrep_x = cfg.nxyz[0]
    expanded: list[list[list[int | float]]] = [[] for _ in range(sim.nocellst)]
    for iz in range(nz):
        base_sim = nx * (nz - iz - 1)
        base_rep = nx * ny * (nz - iz - 1)
        maps = [
            [
                [tgt + (tgt // nrep_x) * nrep_x * (cfg.nxyz[1] - 1), w * weights[0][0]]
                for tgt, w in row
            ]
            for row in cell_ind[base_sim : base_sim + nx]
        ]
        expanded[base_rep : base_rep + nx] = maps
        for j, iy in enumerate(indy[1:]):
            start = base_rep + nx * (j + 1)
            wy = weights[j + 1][0]
            expanded[start : start + nx] = [
                [[tgt + iy * nrep_x, w * wy] for tgt, w in row] for row in maps
            ]
    return expanded


def handle_yaxis_mapping_intensive(
    cfg: BenchmarkConfig,
    sim: SimulationData,
    cell_cent: NDArray,
    refycent: NDArray,
    simycent: NDArray,
) -> NDArray:
    """Extend representative cell indices for y direction (intensive)"""
    indy = np.array([np.argmin(np.abs(simycent - y)) for y in refycent])
    expanded = np.zeros(cfg.nocellsrepgrid, dtype=int)
    nx, ny, nz = cfg.nxyz
    gx = sim.simdim[0]
    for iz in range(nz):
        base_xz = nx * (nz - iz - 1)
        base_xyz = nx * ny * (nz - iz - 1)
        row = cell_cent[base_xz : base_xz + nx]
        mults = np.floor(row / gx) if iz != 0 else 0
        values = row + mults * gx * (sim.simdim[1] - 1)
        expanded[base_xyz : base_xyz + nx] = values
        for j, iy in enumerate(indy[1:]):
            start = base_xyz + nx * (j + 1)
            expanded[start : start + nx] = iy * gx + values
    return expanded


def handle_inactive_mapping(
    cfg: BenchmarkConfig, sim: SimulationData, cell_ind: list[list[list[int | float]]]
) -> NDArray:
    """Identify inactive reporting-grid cells"""
    actindr = []
    for i in sim.actind:
        for mask in cell_ind[i]:
            actindr.append(mask[0])
    actindr = list(dict.fromkeys(actindr))
    allc = np.linspace(0, cfg.nocellsrepgrid - 1, cfg.nocellsrepgrid, dtype=int)
    return np.delete(allc, actindr)


def handle_performance_spatial(
    cfg: BenchmarkConfig,
    sim: SimulationData,
    rstno: list,
    refgrid: tuple[NDArray, NDArray, NDArray, NDArray, NDArray, NDArray],
    mapping: tuple[list, NDArray],
    actindr: NDArray,
) -> None:
    """Create performance spatial maps"""
    cell_ind, _ = mapping
    _, _, _, refxcent, refycent, refzcent = refgrid
    counter = 0.0 * np.ones(cfg.nocellsrepgrid)
    pore_volume = 0.0 * np.ones(cfg.nocellsrepgrid)
    if actindr.size > 0:
        pore_volume[actindr] = 1.0
    latest_dts, cvol_refg, arat_refg, valid = static_map_performance_spatial(
        cfg, sim, cell_ind, counter, pore_volume
    )
    names = ("co2mn", "h2omn", "co2mb", "h2omb")
    for i, rst in enumerate(rstno):
        if i + 1 < len(rstno):
            print(
                f"Processing performance spatial {i+1} out of {len(rstno)}",
                end="\r",
            )
        else:
            print(f"Processing performance spatial {i+1} out of {len(rstno)}")
        arrays = init_performance_arrays(sim, names)
        step_index = rst + sim.noskiprst
        if step_index > 0:
            populate_performance_arrays(sim, arrays, step_index - 1)
        refg = map_performance_to_report_grid(
            cfg, sim, arrays, cell_ind, latest_dts[i], pore_volume, valid
        )
        write_dense_performance_spatial(
            cfg, refg, cvol_refg, arat_refg, refxcent, refycent, refzcent, i
        )


def map_performance_to_report_grid(
    cfg: BenchmarkConfig,
    sim: SimulationData,
    arrays: dict,
    cell_ind: list,
    delta_t: float,
    pore_volume: NDArray,
    valid: NDArray,
) -> dict:
    """Map simulation grid to reporting grid"""
    refg = {
        "co2mn": np.full(cfg.nocellsrepgrid, -np.inf),
        "h2omn": np.full(cfg.nocellsrepgrid, -np.inf),
        "co2mb": np.zeros(cfg.nocellsrepgrid),
        "h2omb": np.zeros(cfg.nocellsrepgrid),
    }
    co2mn = arrays["co2mn_array"]
    h2omn = arrays["h2omn_array"]
    co2mb = arrays["co2mb_array"]
    h2omb = arrays["h2omb_array"]
    ref_co2mn = refg["co2mn"]
    ref_h2omn = refg["h2omn"]
    ref_co2mb = refg["co2mb"]
    ref_h2omb = refg["h2omb"]
    for cell in sim.actind:
        v_co2mn = co2mn[cell]
        v_h2omn = h2omn[cell]
        v_co2mb = co2mb[cell]
        v_h2omb = h2omb[cell]
        for tgt, w in cell_ind[cell]:
            ref_co2mn[tgt] = max(ref_co2mn[tgt], v_co2mn * w)
            ref_h2omn[tgt] = max(ref_h2omn[tgt], v_h2omn * w)
            ref_co2mb[tgt] += v_co2mb * w
            ref_h2omb[tgt] += v_h2omb * w
    ref_co2mn[np.isfinite(ref_co2mn)] *= delta_t
    ref_h2omn[np.isfinite(ref_h2omn)] *= delta_t
    ref_co2mb[valid] = delta_t * ref_co2mb[valid] / pore_volume[valid]
    ref_h2omb[valid] = delta_t * ref_h2omb[valid] / pore_volume[valid]
    return refg


def static_map_performance_spatial(
    cfg: BenchmarkConfig,
    sim: SimulationData,
    cell_ind: list,
    counter: NDArray,
    pore_volume: NDArray,
) -> tuple[list, NDArray, NDArray, NDArray]:
    """Map static quantities for performance spatial data"""
    tmp1, tmp2 = [], []
    with open(
        f"{cfg.flowfol}/{cfg.outfol.split('/')[-1].upper()}.INFOSTEP",
        "r",
        encoding="utf8",
    ) as file:
        for i, row in enumerate(csv.reader(file)):
            if i > 0:
                vals = row[0].strip().split()
                tmp1.append(86400.0 * float(vals[0]))
                tmp2.append(86400.0 * float(vals[1]))
    infotimes = np.array(tmp1)
    tsteps = np.array(tmp2)
    latest_dts = []
    for time_val in cfg.denset[:-1]:
        pos = np.argmin(np.abs(infotimes - (time_val + sim.timeini)))
        latest_dts.append(tsteps[pos - 1] if pos > 0 else 0.0)
    latest_dts.append(tsteps[-1])
    cvol_array = np.zeros(sim.nocellst)
    arat_array = np.zeros(sim.nocellst)
    cvol_refg = np.zeros(cfg.nocellsrepgrid)
    arat_refg = np.zeros(cfg.nocellsrepgrid)
    act = sim.actind
    cvol_array[act] = sim.porva / np.array(sim.init["PORO"])
    if cfg.case != "spe11c":
        arat_array[act] = np.array(sim.init["DZ"]) / np.array(sim.init["DX"])
    else:
        dx = np.array(sim.init["DX"])
        dy = np.array(sim.init["DY"])
        arat_array[act] = np.array(sim.init["DZ"]) / np.sqrt(dx**2 + dy**2)
    for cell in act:
        pv = sim.porv[cell]
        for tgt, _ in cell_ind[cell]:
            cvol_refg[tgt] += cvol_array[cell]
            arat_refg[tgt] += arat_array[cell]
            counter[tgt] += 1.0
            pore_volume[tgt] += pv
    mask = counter > 0.0
    cvol_refg[mask] /= counter[mask]
    arat_refg[mask] /= counter[mask]
    cvol_refg[cvol_refg < 1e-12] = np.nan
    arat_refg[arat_refg < 1e-12] = np.nan
    valid = pore_volume > 0.0
    return latest_dts, cvol_refg, arat_refg, valid


def init_performance_arrays(
    sim: SimulationData, names: tuple[str, str, str, str]
) -> dict:
    """Initialize performance arrays"""
    arrays = {}
    for name in names:
        arrays[f"{name}_array"] = np.zeros(sim.nocellst)
    return arrays


def populate_performance_arrays(
    sim: SimulationData, arrays: dict, step_index: int
) -> None:
    """Populate performance arrays"""
    act = sim.actind
    co2mb = arrays["co2mb_array"]
    h2omb = arrays["h2omb_array"]
    co2mn = arrays["co2mn_array"]
    h2omn = arrays["h2omn_array"]
    porva = sim.porva
    co2mb[act] = np.array(sim.unrst["RES_GAS", step_index + 1])
    if sim.unrst.count("RES_WAT", step_index + 1):
        h2omb[act] = np.array(sim.unrst["RES_WAT", step_index + 1])
    else:
        h2omb[act] = np.array(sim.unrst["RES_OIL", step_index + 1])
    co2mn[act] = np.abs(co2mb[act]) / porva
    h2omn[act] = np.abs(h2omb[act]) / porva


def write_dense_performance_spatial(
    cfg: BenchmarkConfig,
    refg: dict,
    cvol_refg: NDArray,
    arat_refg: NDArray,
    refxcent: NDArray,
    refycent: NDArray,
    refzcent: NDArray,
    i: int,
) -> None:
    """Write dense performance spatial CSV"""
    if cfg.case == "spe11a":
        name_t = f"{round(cfg.denset[i]/3600)}h"
    else:
        name_t = f"{round(cfg.denset[i]/SECONDS_IN_YEAR)}y"
    nx, ny, nz = cfg.nxyz
    co2mn = refg["co2mn"]
    h2omn = refg["h2omn"]
    co2mb = refg["co2mb"]
    h2omb = refg["h2omb"]
    path = f"{cfg.where}/{cfg.case}_performance_spatial_map_{name_t}.csv"
    with open(path, "w", encoding="utf8") as file:
        if cfg.case != "spe11c":
            file.write(
                "# x [m], z [m], cvol [m^3], arat [-], CO2 max_norm_res [-], "
                "H2O max_norm_res [-], CO2 mb_error [-], H2O mb_error [-], post_est [-]"
            )
        else:
            file.write(
                "# x [m], y [m], z [m], cvol [m^3], arat [-], CO2 max_norm_res [-], "
                "H2O max_norm_res [-], CO2 mb_error [-], H2O mb_error [-], post_est [-]"
            )
        idz = 0
        for zcord in refzcent:
            idxy = 0
            basez = -nx * ny * (nz - idz)
            for ycord in refycent:
                for xcord in refxcent:
                    idc = basez + idxy
                    if np.isnan(cvol_refg[idc]):
                        if cfg.case != "spe11c":
                            file.write(
                                f"\n{xcord:.3e}, {zcord:.3e}, n/a, n/a, n/a, n/a, n/a, n/a, n/a"
                            )
                        else:
                            file.write(
                                f"\n{xcord:.3e}, {ycord:.3e}, {zcord:.3e}, n/a, n/a, n/a, n/a, "
                                "n/a, n/a, n/a"
                            )
                    else:
                        if cfg.case != "spe11c":
                            file.write(
                                f"\n{xcord:.3e}, {zcord:.3e}, {cvol_refg[idc]:.3e}, "
                                f"{arat_refg[idc]:.3e}, {co2mn[idc]:.3e}, {h2omn[idc]:.3e}, "
                                f"{co2mb[idc]:.3e}, {h2omb[idc]:.3e}, n/a"
                            )
                        else:
                            file.write(
                                f"\n{xcord:.3e}, {ycord:.3e}, {zcord:.3e}, "
                                f"{cvol_refg[idc]:.3e}, {arat_refg[idc]:.3e}, "
                                f"{co2mn[idc]:.3e}, {h2omn[idc]:.3e}, "
                                f"{co2mb[idc]:.3e}, {h2omb[idc]:.3e}, n/a"
                            )
                    idxy += 1
            idz += 1


def generate_arrays(
    cfg: BenchmarkConfig,
    sim: SimulationData,
    names: list,
    restart_index: int,
    actindr: NDArray,
) -> dict:
    """Populate dense arrays for one restart"""
    arrays = {}
    act = sim.actind
    porva = sim.porva
    for name in names[:-1]:
        arr = np.zeros(sim.nocellst)
        if cfg.case == "spe11a" or (cfg.lower and not sim.cornpoint):
            arr[:] = np.nan
        arrays[f"{name}_array"] = arr
        refg = np.zeros(cfg.nocellsrepgrid)
        if actindr.size > 0:
            refg[actindr] = np.nan
        arrays[f"{name}_refg"] = refg
    tco2_array = np.zeros(sim.nocellst)
    tco2_refg = np.zeros(cfg.nocellsrepgrid)
    sgas = np.abs(np.array(sim.unrst["SGAS", restart_index]))
    rhog = np.array(sim.unrst["GAS_DEN", restart_index])
    pres = np.array(sim.unrst["PRESSURE", restart_index]) - np.array(
        sim.unrst["PCGW", restart_index]
    )
    rhow = np.array(sim.unrst["WAT_DEN", restart_index])
    mask_g = sgas > SGAS_THR
    if not sim.immiscible:
        rss = np.array(sim.unrst["RSW", restart_index])
        rvv = (
            np.array(sim.unrst["RVW", restart_index])
            if sim.unrst.count("RVW", restart_index)
            else 0.0 * rss
        )
        if not sim.isothermal:
            arrays["temp_array"][act] = np.array(sim.unrst["TEMP", restart_index])
        xco2 = rss / (rss + WAT_DEN_REF / GAS_DEN_REF)
        xh2o = rvv / (rvv + GAS_DEN_REF / WAT_DEN_REF)
        co2_g = (1 - xh2o) * sgas * rhog * porva
        co2_d = xco2 * (1 - sgas) * rhow * porva
    else:
        xco2 = 0.0 * sgas
        xh2o = 0.0 * sgas
        co2_g = sgas * rhog * porva
        co2_d = 0.0
    arrays["pressure_array"][act] = 1e5 * pres
    arrays["sgas_array"][act] = sgas * mask_g
    arrays["gden_array"][act] = rhog * mask_g
    arrays["wden_array"][act] = rhow
    arrays["xco2_array"][act] = xco2
    arrays["xh20_array"][act] = xh2o * mask_g
    tco2_array[act] = co2_d + co2_g
    arrays["tco2_array"] = tco2_array
    arrays["tco2_refg"] = tco2_refg
    if cfg.lower and sim.cornpoint:
        pad = np.full(sim.simdim[0] * sim.simdim[1], np.nan)
        for key in arrays:
            if key.endswith("_array") and key != "tco2_array":
                arrays[key] = np.insert(arrays[key], 0, pad)
    return arrays


def map_to_report_grid(
    sim: SimulationData,
    mapping: tuple[list[list[list[int | float]]], NDArray],
    arrays: dict,
) -> None:
    """Map simulation arrays to reporting grid"""
    cell_ind, cell_cent = mapping
    tco2_arr = arrays["tco2_array"]
    tco2_refg = arrays["tco2_refg"]
    for cell in sim.actind:
        tval = tco2_arr[cell]
        for tgt, w in cell_ind[cell]:
            tco2_refg[tgt] += tval * w
    for rep, cent in enumerate(cell_cent):
        cell = int(cent)
        for key in arrays:
            if key.endswith("_array") and not key.startswith("tco2"):
                arrays[key.replace("_array", "_refg")][rep] = arrays[key][cell]


def get_header(
    cfg: BenchmarkConfig, sim: SimulationData, i: int
) -> tuple[str, list[str]]:
    """Get the CSV header for spe11x"""
    if cfg.case == "spe11a":
        name_t = f"{round(cfg.denset[i]/3600)}h"
        text = [
            "# x [m], z [m], pressure [Pa], gas saturation [-], "
            + "mass fraction of CO2 in liquid [-], mass fraction of H20 in vapor [-], "
            + "phase mass density gas [kg/m3], phase mass density water [kg/m3], "
            + "total mass CO2 [kg]"
            + (", temperature [C]" if not sim.isothermal else "")
        ]
    elif cfg.case == "spe11b":
        name_t = f"{round(cfg.denset[i]/SECONDS_IN_YEAR)}y"
        text = [
            "# x [m], z [m], pressure [Pa], gas saturation [-], "
            + "mass fraction of CO2 in liquid [-], mass fraction of H20 in vapor [-], "
            + "phase mass density gas [kg/m3], phase mass density water [kg/m3], "
            + "total mass CO2 [kg], temperature [C]"
        ]
    else:
        name_t = f"{round(cfg.denset[i]/SECONDS_IN_YEAR)}y"
        text = [
            "# x [m], y [m], z [m], pressure [Pa], gas saturation [-], "
            + "mass fraction of CO2 in liquid [-], mass fraction of H20 in vapor [-], "
            + "phase mass density gas [kg/m3], phase mass density water [kg/m3], "
            + "total mass CO2 [kg], temperature [C]"
        ]
    return name_t, text


if __name__ == "__main__":
    main(sys.argv[1:])
