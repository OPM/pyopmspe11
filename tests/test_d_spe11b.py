# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the spe11b case"""

import pathlib
from shutil import copyfile
from opm.io.ecl import ESmry as OpmSummary
from opm.io.ecl import EclFile as OpmFile
from pyopmspe11.core.pyopmspe11 import main

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_d_spe11b(test_d_f_g_workdir, monkeypatch):
    """Run spe11b via main() and validate outputs"""
    copyfile(testpth / "configs/input.toml", test_d_f_g_workdir / "input.toml")
    monkeypatch.chdir(test_d_f_g_workdir)
    main()
    case = test_d_f_g_workdir / "output/flow/OUTPUT"
    assert case.with_suffix(".UNRST").exists()
    assert abs(OpmSummary(f"{case}.SMSPEC")["FGMIP"][-1] - 8.2760104e07) < 1e5
    assert abs(sum(OpmFile(f"{case}.INIT")["PORV"]) - 2.0512318e07) < 1e-6
