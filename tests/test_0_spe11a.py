# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the spe11a case"""

import os
import pathlib
import subprocess

dirname: pathlib.Path = pathlib.Path(__file__).parent


def test_spe11a():
    """See configs/spe11a_data_format.txt"""
    os.chdir(f"{dirname}/configs")
    subprocess.run(
        [
            "pyopmspe11",
            "-i",
            "spe11a_data_format.txt",
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
            "-u",
            "opm",
        ],
        check=True,
    )
    assert os.path.exists(f"{dirname}/configs/spe11a/data/spe11a_time_series.csv")
