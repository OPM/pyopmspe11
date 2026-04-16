# SPDX-FileCopyrightText: 2024-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test keywords in the main deck, including source positions"""

import pathlib
import subprocess
import filecmp

testpth: pathlib.Path = pathlib.Path(__file__).parent
mainpth: pathlib.Path = pathlib.Path(__file__).parents[1]


def test_b_deck_content(tmp_path, monkeypatch):
    """Compare generated main decks against references"""
    monkeypatch.chdir(tmp_path)
    for x, y, size in zip(
        ("A", "B", "C"), ("a", "b", "c"), ("1cm", "10m", "50m-50m-10m")
    ):
        subprocess.run(
            [
                "pyopmspe11",
                "-i",
                f"{mainpth}/benchmark/spe11{y}/r1_Cart_{size}.toml",
                "-o",
                f"spe11{y}_deck",
                "-m",
                "deck",
                "-f",
                "0",
            ],
            check=True,
        )
        ref = testpth / f"decks/SPE11{x}.DATA"
        out = tmp_path / f"spe11{y}_deck/SPE11{x}_DECK.DATA"
        assert filecmp.cmp(
            ref, out, shallow=False
        ), f"Difference between SPE11{x}_DECK.DATA and generated deck"
