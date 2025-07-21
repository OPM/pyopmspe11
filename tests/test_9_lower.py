# SPDX-FileCopyrightText: 2025 NORCE
# SPDX-License-Identifier: MIT

"""Test the generation of the lower localized domain (flag '-n lower')"""

import pathlib
import subprocess
import numpy as np

testpth: pathlib.Path = pathlib.Path(__file__).parent
mainpth: pathlib.Path = pathlib.Path(__file__).parents[1]


def test_lower():
    """See examples/spe11x.toml"""
    cases = ["a", "b", "c"]
    reportings = ["15,1,35", "15,1,35", "15,5,35"]
    times = ["h", "y", "y"]
    nans = [358, 327, 1460]
    for c, r, t, n in zip(cases, reportings, times, nans):
        subprocess.run(
            [
                "pyopmspe11",
                "-i",
                f"{mainpth}/examples/spe11{c}.toml",
                "-o",
                f"{testpth}/decks/lower_spe11{c}",
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
            f"{testpth}/decks/lower_spe11{c}/spe11{c}_spatial_map_0{t}.csv",
            delimiter=",",
            skip_header=1,
        )
        assert np.sum(np.isnan(csv[:, 3])) == n
