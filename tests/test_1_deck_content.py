# SPDX-FileCopyrightText: 2024 NORCE
# SPDX-License-Identifier: MIT

"""Test the keywords in the main deck, including source positions"""

import filecmp
import pathlib
import subprocess

testpth: pathlib.Path = pathlib.Path(__file__).parent
mainpth: pathlib.Path = pathlib.Path(__file__).parents[1]


def test_deck_content():
    """See benchmark/spe11x/r1_x (cases with reporting grids)"""
    for x, size in zip(["A", "B", "C"], ["1cm", "10m", "50m-50m-10m"]):
        subprocess.run(
            [
                "pyopmspe11",
                "-i",
                f"{mainpth}/benchmark/spe11{x.lower()}/r1_Cart_{size}.txt",
                "-o",
                f"{testpth}/decks/SPE11{x}",
                "-m",
                "deck",
            ],
            check=True,
        )
        text = f"Difference between SPE11{x}.DATA and {testpth}/decks/SPE11{x}/deck/SPE11{x}.DATA"
        assert filecmp.cmp(
            f"{testpth}/decks/SPE11{x}.DATA",
            f"{testpth}/decks/SPE11{x}/deck/SPE11{x}.DATA",
            shallow=False,
        ), text
