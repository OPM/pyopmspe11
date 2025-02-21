# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the spe11c case"""

import os
import pathlib
import subprocess

dirname: pathlib.Path = pathlib.Path(__file__).parent


def test_spe11c():
    """See configs/spe11c.txt"""
    os.chdir(f"{dirname}/configs")
    subprocess.run(
        [
            "pyopmspe11",
            "-o",
            "spe11c",
            "-i",
            "spe11c.txt",
            "-m",
            "all",
            "-g",
            "dense_performance",
            "-r",
            "24,3,12",
            "-t",
            "5",
        ],
        check=True,
    )
    assert os.path.exists(f"{dirname}/configs/spe11c/figures/spe11c_temp_2Dmaps.png")
