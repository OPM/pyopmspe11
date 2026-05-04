# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the spe11c case"""

import pathlib
import subprocess

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_e_spe11c(tmp_path, monkeypatch):
    """Run spe11c and validate outputs"""
    monkeypatch.chdir(tmp_path)
    spe11c = (testpth / "configs/spe11c_data_format.toml").resolve()
    subprocess.run(
        [
            "pyopmspe11",
            "-m",
            "deck_flow_data",
            "-g",
            "dense_performance_sparse",
            "-w",
            "0.1",
            "-i",
            str(spe11c),
            "-o",
            "spe11c",
            "-r",
            "168,100,120",
            "-t",
            "0,5,10,15,20,25,30,35,40,45,50,75,100,150,200,250,300,350,400,450,500,"
            + "600,700,800,900,1000",
            "-w",
            "0.1",
        ],
        check=True,
    )
    for file in ("check_format", "is_notebook"):
        subprocess.run(
            [
                "curl",
                "-o",
                f"{file}.py",
                "https://raw.githubusercontent.com/Simulation-Benchmarks/11thSPE-CSP/"
                + f"main/evaluation/{file}.py",
            ],
            check=True,
        )
    result = subprocess.run(
        ["python3", "check_format.py", "-f", "./spe11c/data", "-c", "C"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert result.stdout.count("Successfully") == 2, "Issue in spe11c/data"
