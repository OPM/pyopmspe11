# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the scrip to compare the data for different runs"""

import os


def test_comparison():
    """See visualization/plot.py"""
    cwd = os.getcwd()
    os.system(f"mkdir {cwd}/tests/configs/test_compare")
    os.system(f"mkdir {cwd}/tests/configs/test_compare/spe11c_ens1")
    os.system(f"mkdir {cwd}/tests/configs/test_compare/spe11c_ens2")
    os.system(
        f"cp -R {cwd}/tests/configs/spe11c/. {cwd}/tests/configs/test_compare/spe11c_ens1"
    )
    os.system(
        f"cp -R {cwd}/tests/configs/spe11c/. {cwd}/tests/configs/test_compare/spe11c_ens2"
    )
    os.chdir(f"{cwd}/tests/configs/test_compare")
    os.system("pyopmspe11 -c spe11c -m performance")
    os.chdir(f"{cwd}/tests/configs")
    assert os.path.exists(
        f"{cwd}/tests/configs/test_compare/compare/spe11c_performance_detailed.png"
    )
    assert os.path.exists(
        f"{cwd}/tests/configs/test_compare/compare/spe11c_performance.png"
    )
    os.chdir(cwd)
