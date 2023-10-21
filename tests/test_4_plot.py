# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the scrip to plot the data"""

import os
from pyopmspe11.visualization.plotting import main


def test_plot():
    """See visualization/plotting.py"""
    cwd = os.getcwd()
    os.chdir(f"{cwd}/tests/configs")
    os.mkdir(f"{cwd}/tests/configs/output/figures")
    main()
    assert os.path.exists(f"{cwd}/tests/configs/output/figures/spe11b_sparse_data.png")
    os.chdir(cwd)
