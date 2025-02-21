# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the scrip to plot the data"""

import os
import pathlib
from pyopmspe11.visualization.plotting import main

dirname: pathlib.Path = pathlib.Path(__file__).parent


def test_plot():
    """See visualization/plotting.py"""
    message = "Please run first test_1_spe11b"
    assert os.path.exists(f"{dirname}/configs/output"), message
    message = "Please run first test_3_data"
    assert os.path.exists(f"{dirname}/configs/output/data"), message
    os.chdir(f"{dirname}/configs")
    os.mkdir(f"{dirname}/configs/output/figures")
    main()
    assert os.path.exists(f"{dirname}/configs/output/figures/spe11b_sparse_data.png")
