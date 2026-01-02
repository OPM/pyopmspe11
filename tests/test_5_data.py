# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the scrip to write the data as required in the benchmark"""

import os
import pathlib
from pyopmspe11.visualization.data import main

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_data():
    """See src/pyopmspe11/visualization/data.py"""
    message = "Please run first test_3_spe11b"
    assert os.path.exists(f"{testpth}/output/output"), message
    if not os.path.exists(f"{testpth}/output/output/data"):
        os.mkdir(f"{testpth}/output/output/data")
    os.chdir(f"{testpth}/output")
    main()
    assert os.path.exists(f"{testpth}/output/output/data/spe11b_time_series.csv")
