# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the spe11c case"""

import pathlib
import subprocess
from opm.io.ecl import ESmry as OpmSummary
from opm.io.ecl import EclFile as OpmFile

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_e_spe11c(test_e_h_workdir, monkeypatch):
    """Run spe11c and validate outputs"""
    monkeypatch.chdir(test_e_h_workdir)
    subprocess.run(
        [
            "pyopmspe11",
            "-o",
            "spe11c",
            "-i",
            f"{testpth}/configs/spe11c.toml",
            "-m",
            "all",
            "-g",
            "dense_performance_sparse",
            "-r",
            "24,3,12",
            "-t",
            "25",
        ],
        check=True,
    )
    fig = test_e_h_workdir / "spe11c/figures/spe11c_temp_2Dmaps.png"
    assert fig.exists()
    case = test_e_h_workdir / "spe11c/flow/SPE11C"
    assert abs(OpmSummary(f"{case}.SMSPEC")["FGMIP"][-1] - 1.1813521e11) < 1e8
    assert abs(sum(OpmFile(f"{case}.INIT")["PORV"]) - 2.6425968e11) < 1e-6
