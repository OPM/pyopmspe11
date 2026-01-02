# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the spe11a case"""

import os
import pathlib
import subprocess
from opm.io.ecl import ESmry as OpmSummary
from opm.io.ecl import EclFile as OpmFile

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_spe11a():
    """See configs/spe11a_data_format.toml"""
    if not os.path.exists(f"{testpth}/output"):
        os.mkdir(f"{testpth}/output")
    if not os.path.exists(f"{testpth}/output/benchmark"):
        os.mkdir(f"{testpth}/output/benchmark")
    os.chdir(f"{testpth}/output/benchmark")
    subprocess.run(
        [
            "pyopmspe11",
            "-i",
            f"{testpth}/configs/spe11a_data_format.toml",
            "-o",
            "spe11a",
            "-m",
            "deck_flow_data",
            "-g",
            "all",
            "-r",
            "280,1,120",
            "-t",
            "1",
            "-w",
            "0.16666666666666666",
        ],
        check=True,
    )
    file = f"{testpth}/output/benchmark/spe11a/data/spe11a_time_series.csv"
    assert os.path.exists(file)
    case = f"{testpth}/output/benchmark/spe11a/flow/SPE11A"
    assert abs(OpmSummary(f"{case}.SMSPEC")["FGMIP"][-1] - 0.0046843435) < 1e-6
    assert abs(sum(OpmFile(f"{case}.INIT")["PORV"]) - 0.013631029) < 1e-6
