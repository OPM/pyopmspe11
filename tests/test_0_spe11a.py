# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the spe11a case"""

import os
import subprocess
from subprocess import PIPE, Popen


def test_spe11a():
    """See configs/spe11a.txt"""
    cwd = os.getcwd()
    os.chdir(f"{cwd}/tests/configs")
    subprocess.run(
        [
            "pyopmspe11",
            "-i",
            "spe11a.txt",
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
    subprocess.run(
        [
            "curl",
            "-o",
            "./check_format.py",
            "https://raw.githubusercontent.com/Simulation-Benchmarks/11thSPE-CSP/"
            + "main/evaluation/check_format.py",
        ],
        check=True,
    )
    with Popen(
        args="python3 check_format.py -f ./spe11a/data -c A", stdout=PIPE, shell=True
    ) as process:
        check = str(process.communicate()[0])[4:-3]
    assert check.count("Successfully") == 2
    os.chdir(cwd)
