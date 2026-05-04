# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the script to write the data as required in the benchmark.
Note we use smaller grids and simulation times for testing purposes"""

import pathlib
import filecmp
import shutil

from pyopmspe11.visualization.data import main

testpth = pathlib.Path(__file__).parent
mainpth = pathlib.Path(__file__).parents[1]


def test_f_data(tmp_path, monkeypatch):
    """Generate benchmark data output"""
    monkeypatch.chdir(tmp_path)
    cases = {
        "a": {"t": "1", "w": "0.16666666666666666", "r": "28,1,12"},
        "b": {"t": "5", "w": "0.1", "r": "84,1,12"},
        "c": {"t": "5", "w": "0.1", "r": "16,10,12"},
    }
    for x, flagspec in cases.items():
        for grid in ["cartesian", "corner-point"]:
            run = tmp_path / f"spe11{x}_{grid}"
            shutil.copytree(testpth / "flows" / f"spe11{x}_{grid}", run)
            flags = [
                "-p",
                f"spe11{x}_{grid}",
                "-f",
                "0",
                "-g",
                "all",
                "-d",
                f"spe11{x}",
                "-w",
                flagspec["w"],
                "-r",
                flagspec["r"],
                "-t",
                flagspec["t"],
            ]
            main(flags)
            ref = testpth / "datas" / f"spe11{x}_{grid}"
            files = sorted(p.name for p in ref.iterdir() if p.is_file())
            match, mismatch, error = filecmp.cmpfiles(run, ref, files, shallow=False)
            assert (
                match == files
            ), f"Config mismatch for spe11{x}_{grid}: mismatch={mismatch}, error={error}"
