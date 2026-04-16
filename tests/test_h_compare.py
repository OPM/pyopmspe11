# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the script to compare the data for different runs"""

import subprocess
from shutil import copytree


def test_h_compare(test_e_h_workdir, monkeypatch):
    """Compare performance plots between two spe11c runs"""
    src = test_e_h_workdir / "spe11c"
    assert src.exists(), "Please run test_e_spe11c first"
    testcmp = test_e_h_workdir / "test_compare"
    ens1 = testcmp / "spe11c_ens1"
    ens2 = testcmp / "spe11c_ens2"
    copytree(src, ens1)
    copytree(src, ens2)
    monkeypatch.chdir(testcmp)
    subprocess.run(["pyopmspe11", "-c", "spe11c", "-m", "performance"], check=True)
    cmpdir = testcmp / "compare"
    assert (cmpdir / "spe11c_performance_detailed.png").exists()
    assert (cmpdir / "spe11c_performance.png").exists()
