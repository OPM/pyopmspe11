# SPDX-FileCopyrightText: 2024-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the format of the generated data for the three cases"""

import pathlib
import subprocess

testpth: pathlib.Path = pathlib.Path(__file__).parent.resolve()


def test_i_data_format(test_c_i_workdir, monkeypatch):
    """Run official spe11 format checks"""
    run = test_c_i_workdir
    monkeypatch.chdir(run)
    assert (run / "spe11a").exists(), "Please run test_c_spe11a first"
    spe11b = (testpth / "configs/spe11b_data_format.toml").resolve()
    spe11c = (testpth / "configs/spe11c_data_format.toml").resolve()
    base = [
        "pyopmspe11",
        "-m",
        "deck_flow_data",
        "-g",
        "dense_performance_sparse",
        "-w",
        "0.1",
    ]
    subprocess.run(
        base + ["-i", str(spe11b), "-o", "spe11b", "-r", "840,1,120", "-t", "5"],
        check=True,
    )
    subprocess.run(
        base
        + [
            "-i",
            str(spe11c),
            "-o",
            "spe11c",
            "-r",
            "168,100,120",
            "-t",
            "0,5,10,15,20,25,30,35,40,45,50,75,100,150,200,250,300,350,400,450,500,"
            + "600,700,800,900,1000",
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
    for x in ("a", "b", "c"):
        result = subprocess.run(
            ["python3", "check_format.py", "-f", f"./spe11{x}/data", "-c", x.upper()],
            capture_output=True,
            text=True,
            check=True,
        )
        assert result.stdout.count("Successfully") == 2, f"Issue in spe11{x}/data"
