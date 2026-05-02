# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the script to compare the data for different runs"""

import pathlib
from shutil import copytree

from pyopmspe11.core.pyopmspe11 import main

testpth = pathlib.Path(__file__).parent


def test_h_compare(tmp_path, monkeypatch):
    """Run compare via main()"""
    monkeypatch.chdir(tmp_path)
    for x in ("cartesian", "corner-point"):
        run = tmp_path / f"spe11c_{x}"
        copytree(testpth / "datas" / f"spe11c_{x}", run)
    main(["-c", "spe11c"])
    for file in ["performance", "performance_detailed", "sparse_data"]:
        figure = tmp_path / "compare" / f"spe11c_{file}.png"
        assert figure.is_file(), f"Missing {figure}"
