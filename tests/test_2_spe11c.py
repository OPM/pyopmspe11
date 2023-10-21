# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the spe11c case"""

import os
import subprocess


def test_spe11b():
    """See configs/spe11c.txt"""
    cwd = os.getcwd()
    os.chdir(f"{os.getcwd()}/tests/configs")
    subprocess.run(
        [
            "pyopmspe11",
            "-i",
            "spe11c.txt",
            "-o",
            "spe11c",
            "-m",
            "all",
            "-g",
            "all",
            "-r",
            "50,5,25",
            "-t",
            "5",
        ],
        check=True,
    )
    assert os.path.exists(f"{cwd}/tests/configs/spe11c/figures/spe11c_temp_2Dmaps.png")
    os.chdir(cwd)
