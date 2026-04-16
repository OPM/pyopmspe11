# SPDX-FileCopyrightText: 2024-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test parsing equivalence of .txt and .toml configuration files"""

import filecmp
import pathlib
import subprocess

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_a_txt_toml(tmp_path, monkeypatch):
    """Compare generated decks from .txt and .toml configs"""
    monkeypatch.chdir(tmp_path)
    for spe in ("spe11a", "spe11b", "spe11c"):
        decks = []
        for ext in ("txt", "toml"):
            subprocess.run(
                [
                    "pyopmspe11",
                    "-i",
                    f"{testpth}/configs/{spe}_data_format.{ext}",
                    "-o",
                    f"{spe}_{ext}",
                    "-m",
                    "deck",
                ],
                check=True,
            )
            src = tmp_path / f"{spe}_{ext}/deck/{spe.upper()}_{ext.upper()}.DATA"
            dst = tmp_path / f"{spe}_{ext}/deck/{spe.upper()}.DATA"
            src.replace(dst)
            decks.append(tmp_path / f"{spe}_{ext}/deck")
        files = [f"{spe.upper()}.DATA", "TABLES.INC", "FIPNUM.INC"]
        if spe != "spe11a":
            files.append("PVBOUNDARIES.INC")
        if spe == "spe11c":
            files.append("GRID.INC")
        match, mismatch, error = filecmp.cmpfiles(
            decks[0], decks[1], files, shallow=False
        )
        assert (
            match == files
        ), f"Config mismatch for {spe}: mismatch={mismatch}, error={error}"
