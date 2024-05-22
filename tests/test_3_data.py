# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the scrip to write the data as required in the benchmark"""

import os
from pyopmspe11.visualization.data import main


def test_data():
    """See visualization/data.py"""
    cwd = os.getcwd()
    os.mkdir(f"{cwd}/tests/configs/output/data")
    os.chdir(f"{cwd}/tests/configs")
    main()
    assert os.path.exists(f"{cwd}/tests/configs/output/data/spe11b_time_series.csv")
    os.chdir(cwd)
