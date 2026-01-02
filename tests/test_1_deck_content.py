# SPDX-FileCopyrightText: 2024-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the keywords in the main deck, including source positions"""

import os
import filecmp
import pathlib
import subprocess

testpth: pathlib.Path = pathlib.Path(__file__).parent
mainpth: pathlib.Path = pathlib.Path(__file__).parents[1]


def test_deck_content():
    """See benchmark/spe11x/r1_x (cases with reporting grids)"""
    if not os.path.exists(f"{testpth}/output"):
        os.mkdir(f"{testpth}/output")
    for x, y, size in zip(
        ["A", "B", "C"], ["a", "b", "c"], ["1cm", "10m", "50m-50m-10m"]
    ):
        subprocess.run(
            [
                "pyopmspe11",
                "-i",
                f"{mainpth}/benchmark/spe11{y}/r1_Cart_{size}.toml",
                "-o",
                f"{testpth}/output/spe11{y}_deck",
                "-m",
                "deck",
                "-f",
                "0",
            ],
            check=True,
        )
        assert filecmp.cmp(
            f"{testpth}/decks/SPE11{x}.DATA",
            f"{testpth}/output/spe11{y}_deck/SPE11{x}_DECK.DATA",
            shallow=False,
        ), (
            f"Difference between SPE11{x}_DECK.DATA and "
            f"{testpth}/decks/spe11{y}_deck/SPE11{x}_DECK.DATA"
        )
