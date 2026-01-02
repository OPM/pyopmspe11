# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the scrip to compare the data for different runs"""

import os
import pathlib

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_comparison():
    """See src/pyopmspe11/visualization/plotting.py"""
    message = "Please run first test_4_spe11c"
    assert os.path.exists(f"{testpth}/output/spe11c"), message
    os.system(f"mkdir {testpth}/output/test_compare")
    os.system(f"mkdir {testpth}/output/test_compare/spe11c_ens1")
    os.system(f"mkdir {testpth}/output/test_compare/spe11c_ens2")
    os.system(
        f"cp -R {testpth}/output/spe11c/. {testpth}/output/test_compare/spe11c_ens1"
    )
    os.system(
        f"cp -R {testpth}/output/spe11c/. {testpth}/output/test_compare/spe11c_ens2"
    )
    os.chdir(f"{testpth}/output/test_compare")
    os.system("pyopmspe11 -c spe11c -m performance")
    os.chdir(f"{testpth}/output")
    assert os.path.exists(
        f"{testpth}/output/test_compare/compare/spe11c_performance_detailed.png"
    )
    assert os.path.exists(
        f"{testpth}/output/test_compare/compare/spe11c_performance.png"
    )
