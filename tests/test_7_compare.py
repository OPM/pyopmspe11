# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the scrip to compare the data for different runs"""

import os
import pathlib

dirname: pathlib.Path = pathlib.Path(__file__).parent


def test_comparison():
    """See visualization/plot.py"""
    message = "Please run first test_2_spe11c"
    assert os.path.exists(f"{dirname}/configs/spe11c"), message
    os.system(f"mkdir {dirname}/configs/test_compare")
    os.system(f"mkdir {dirname}/configs/test_compare/spe11c_ens1")
    os.system(f"mkdir {dirname}/configs/test_compare/spe11c_ens2")
    os.system(
        f"cp -R {dirname}/configs/spe11c/. {dirname}/configs/test_compare/spe11c_ens1"
    )
    os.system(
        f"cp -R {dirname}/configs/spe11c/. {dirname}/configs/test_compare/spe11c_ens2"
    )
    os.chdir(f"{dirname}/configs/test_compare")
    os.system("pyopmspe11 -c spe11c -m performance")
    os.chdir(f"{dirname}/configs")
    assert os.path.exists(
        f"{dirname}/configs/test_compare/compare/spe11c_performance_detailed.png"
    )
    assert os.path.exists(
        f"{dirname}/configs/test_compare/compare/spe11c_performance.png"
    )
