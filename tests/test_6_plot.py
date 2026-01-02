# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the scrip to plot the data"""

import os
import pathlib
from pyopmspe11.visualization.plotting import main

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_plot():
    """See src/pyopmspe11/visualization/plotting.py"""
    message = "Please run first test_1_spe11b"
    assert os.path.exists(f"{testpth}/output/output"), message
    message = "Please run first test_3_data"
    assert os.path.exists(f"{testpth}/output/output/data"), message
    os.chdir(f"{testpth}/output")
    if not os.path.exists(f"{testpth}/output/output/figures"):
        os.mkdir(f"{testpth}/output/output/figures")
    main()
    assert os.path.exists(f"{testpth}/output/output/figures/spe11b_sparse_data.png")
