# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the spe11a case"""

import pathlib
import subprocess
from opm.io.ecl import ESmry as OpmSummary
from opm.io.ecl import EclFile as OpmFile

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_c_spe11a(test_c_i_workdir, monkeypatch):
    """Run spe11a and validate outputs"""
    monkeypatch.chdir(test_c_i_workdir)
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
    csv = test_c_i_workdir / "spe11a/data/spe11a_time_series.csv"
    assert csv.exists()
    case = test_c_i_workdir / "spe11a/flow/SPE11A"
    assert abs(OpmSummary(f"{case}.SMSPEC")["FGMIP"][-1] - 0.0046843435) < 1e-6
    assert abs(sum(OpmFile(f"{case}.INIT")["PORV"]) - 0.013631029) < 1e-6
