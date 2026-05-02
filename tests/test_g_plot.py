# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the script to plot the data"""

import pathlib
import shutil

from pyopmspe11.visualization.plotting import main

testpth = pathlib.Path(__file__).parent


def test_g_plot(tmp_path, monkeypatch):
    """Generate benchmark plots"""
    monkeypatch.chdir(tmp_path)
    tflags = {"a": "1", "b": "5", "c": "5"}
    for x in ("a", "b", "c"):
        run = tmp_path / f"spe11{x}_corner-point"
        shutil.copytree(testpth / "datas" / f"spe11{x}_corner-point", run)
        flags = [
            "-p",
            f"spe11{x}_corner-point",
            "-g",
            "all",
            "-f",
            "0",
            "-d",
            f"spe11{x}",
            "-t",
            tflags[x],
        ]
        main(flags)
        for file in [
            "performance",
            "performance_detailed",
            "sparse_data",
            "CO2 max_norm_res_2Dmaps",
            "CO2 mb_error_2Dmaps",
            "H2O max_norm_res_2Dmaps",
            "H2O mb_error_2Dmaps",
            "arat_2Dmaps",
            "cvol_2Dmaps",
            "gden_2Dmaps",
            "pressure_2Dmaps",
            "sgas_2Dmaps",
            "tco2_2Dmaps",
            "wden_2Dmaps",
            "xco2_2Dmaps",
            "xh20_2Dmaps",
        ]:
            figure = run / f"spe11{x}_{file}.png"
            assert figure.is_file(), f"Missing {figure}"
