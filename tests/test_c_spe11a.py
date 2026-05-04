# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the spe11a case"""

import pathlib
import subprocess

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_c_spe11a(tmp_path, monkeypatch):
    """Run spe11a and validate outputs"""
    monkeypatch.chdir(tmp_path)
    spe11a = (testpth / "configs/spe11a_data_format.toml").resolve()
    subprocess.run(
        [
            "pyopmspe11",
            "-g",
            "dense_performance_sparse",
            "-m",
            "deck_flow_data",
            "-w",
            "0.1",
            "-i",
            str(spe11a),
            "-o",
            "spe11a",
            "-r",
            "280,1,120",
            "-t",
            "1",
            "-w",
            "0.16666666666666666",
        ],
        check=True,
    )
    for file in ("check_format", "is_notebook"):
        subprocess.run(
            [
                "curl",
                "https://raw.githubusercontent.com/Simulation-Benchmarks/11thSPE-CSP/"
                + f"main/evaluation/{file}.py",
                "-o",
                f"{file}.py",
            ],
            check=True,
        )
    result = subprocess.run(
        ["python3", "check_format.py", "-f", "./spe11a/data", "-c", "A"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert result.stdout.count("Successfully") == 2, "Issue in spe11a/data"
