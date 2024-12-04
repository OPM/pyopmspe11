# SPDX-FileCopyrightText: 2024 NORCE
# SPDX-License-Identifier: MIT

"""Test the format of the generated data for the three cases"""

import os
import pathlib
import subprocess
from subprocess import PIPE, Popen

dirname: pathlib.Path = pathlib.Path(__file__).parent


def test_format():
    """See https://github.com/Simulation-Benchmarks/11thSPE-CSP/blob/
    main/evaluation/check_format.py"""
    message = "Please run first test_0_spe11a"
    assert os.path.exists(f"{dirname}/configs/spe11a"), message
    os.chdir(f"{dirname}/configs")
    if os.path.exists(f"{dirname}/configs/spe11c"):
        os.system(f"rm -rf {dirname}/configs/spe11c")
    subprocess.run(
        [
            "pyopmspe11",
            "-i",
            "spe11b_data_format.txt",
            "-o",
            "spe11b",
            "-m",
            "deck_flow_data",
            "-g",
            "dense_performance_sparse",
            "-r",
            "840,1,120",
            "-t",
            "5",
            "-w",
            "0.1",
            "-u",
            "opm",
        ],
        check=True,
    )
    subprocess.run(
        [
            "pyopmspe11",
            "-i",
            "spe11c_data_format.txt",
            "-o",
            "spe11c",
            "-m",
            "deck_flow_data",
            "-g",
            "dense_performance_sparse",
            "-r",
            "168,100,120",
            "-t",
            "0,5,10,15,20,25,30,35,40,45,50,75,100,150,200,"
            + "250,300,350,400,450,500,600,700,800,900,1000",
            "-w",
            "0.1",
            "-u",
            "opm",
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
    with Popen(
        args="python3 check_format.py -f ./spe11b/data -c B", stdout=PIPE, shell=True
    ) as process:
        check = str(process.communicate()[0])[4:-3]
    assert check.count("Successfully") == 2
    with Popen(
        args="python3 check_format.py -f ./spe11c/data -c C", stdout=PIPE, shell=True
    ) as process:
        check = str(process.communicate()[0])[4:-3]
    assert check.count("Successfully") == 2
