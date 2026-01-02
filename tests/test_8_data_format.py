# SPDX-FileCopyrightText: 2024-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the format of the generated data for the three cases"""

import os
import pathlib
import subprocess
from subprocess import PIPE, Popen

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_format():
    """See https://github.com/Simulation-Benchmarks/11thSPE-CSP/blob/
    main/evaluation/check_format.py"""
    message = "Please run first test_2_spe11a"
    assert os.path.exists(f"{testpth}/output/benchmark/spe11a"), message
    os.chdir(f"{testpth}/output/benchmark")
    spe11b = f"{testpth}/configs/spe11b_data_format.toml"
    spe11c = f"{testpth}/configs/spe11c_data_format.toml"
    base = "pyopmspe11 -m deck_flow_data -g dense_performance_sparse -w 0.1 "
    command = base + f"-i {spe11b} -o spe11b -r 840,1,120 -t 5 & "
    command += base + f"-i {spe11c} -o spe11c -r 168,100,120 -t "
    command += "0,5,10,15,20,25,30,35,40,45,50,75,100,150,200,"
    command += "250,300,350,400,450,500,600,700,800,900,1000 & wait"
    os.system(command)
    for file in ["check_format", "is_notebook"]:
        subprocess.run(
            [
                "curl",
                "-o",
                f"./{file}.py",
                "https://raw.githubusercontent.com/Simulation-Benchmarks/11thSPE-CSP/"
                + f"main/evaluation/{file}.py",
            ],
            check=True,
        )
    for x in ["a", "b", "c"]:
        arg = f"python3 check_format.py -f ./spe11{x}/data -c {x.upper()}"
        with Popen(args=arg, stdout=PIPE, shell=True) as process:
            check = str(process.communicate()[0])[4:-3]
        assert check.count("Successfully") == 2, f"Issue in spe11{x}/data"
