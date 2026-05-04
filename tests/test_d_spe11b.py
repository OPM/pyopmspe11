# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the spe11b case"""

import pathlib
import subprocess

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_d_spe11b(tmp_path, monkeypatch):
    """Run spe11c and validate outputs"""
    monkeypatch.chdir(tmp_path)
    spe11b = (testpth / "configs/spe11b_data_format.toml").resolve()
    subprocess.run(
        [
            "pyopmspe11",
            "-w",
            "0.1",
            "-m",
            "deck_flow_data",
            "-g",
            "dense_performance_sparse",
            "-i",
            str(spe11b),
            "-o",
            "spe11b",
            "-r",
            "840,1,120",
            "-t",
            "5",
            "-w",
            "0.1",
        ],
        check=True,
    )
    for file in ("is_notebook", "check_format"):
        subprocess.run(
            [
                "curl",
                "-o",
                f"{file}.py",
                "https://raw.githubusercontent.com/Simulation-Benchmarks/"
                + f"11thSPE-CSP/main/evaluation/{file}.py",
            ],
            check=True,
        )
    result = subprocess.run(
        ["python3", "check_format.py", "-f", "./spe11b/data", "-c", "B"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert result.stdout.count("Successfully") == 2, "Issue in spe11b/data"
