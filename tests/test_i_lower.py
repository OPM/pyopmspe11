# SPDX-FileCopyrightText: 2025-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test generation of the lower localized domain (flag '-n lower')"""

import pathlib
import subprocess
import numpy as np

testpth: pathlib.Path = pathlib.Path(__file__).parent
mainpth: pathlib.Path = pathlib.Path(__file__).parents[1]


def test_i_lower(tmp_path, monkeypatch):
    """Run lower-domain cases and validate NaN counts"""
    monkeypatch.chdir(tmp_path)
    cases = ("a", "b", "c")
    reportings = ("15,1,35", "15,1,35", "15,5,35")
    times = ("h", "y", "y")
    nans = (358, 327, 1460)
    for c, r, t, n in zip(cases, reportings, times, nans):
        out = tmp_path / f"lower_spe11{c}"
        subprocess.run(
            [
                "pyopmspe11",
                "-i",
                f"{mainpth}/examples/spe11{c}.toml",
                "-o",
                out.name,
                "-m",
                "deck_flow_data",
                "-f",
                "0",
                "-g",
                "dense",
                "-r",
                r,
                "-n",
                "lower",
            ],
            check=True,
        )
        csv = np.genfromtxt(
            out / f"spe11{c}_spatial_map_0{t}.csv", delimiter=",", skip_header=1
        )
        assert np.sum(np.isnan(csv[:, 3])) == n
