# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the spe11b case"""

import os
import pathlib
from opm.io.ecl import ESmry as OpmSummary
from opm.io.ecl import EclFile as OpmFile
from pyopmspe11.core.pyopmspe11 import main

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_spe11b():
    """See configs/input.toml"""
    if not os.path.exists(f"{testpth}/output"):
        os.mkdir(f"{testpth}/output")
    os.chdir(f"{testpth}/output")
    os.system(f"cp {testpth}/configs/input.toml .")
    main()
    assert os.path.exists(f"{testpth}/output/output/flow/OUTPUT.UNRST")
    case = f"{testpth}/output/output/flow/OUTPUT"
    assert abs(OpmSummary(f"{case}.SMSPEC")["FGMIP"][-1] - 8.2760104e07) < 1e5
    assert abs(sum(OpmFile(f"{case}.INIT")["PORV"]) - 2.0512318e07) < 1e-6
