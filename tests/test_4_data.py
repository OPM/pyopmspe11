# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the scrip to write the data as required in the benchmark"""

import os
import pathlib
from pyopmspe11.visualization.data import main

dirname: pathlib.Path = pathlib.Path(__file__).parent


def test_data():
    """See visualization/data.py"""
    message = "Please run first test_1_spe11b"
    assert os.path.exists(f"{dirname}/configs/output"), message
    os.mkdir(f"{dirname}/configs/output/data")
    os.chdir(f"{dirname}/configs")
    main()
    assert os.path.exists(f"{dirname}/configs/output/data/spe11b_time_series.csv")
