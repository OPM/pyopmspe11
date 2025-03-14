# SPDX-FileCopyrightText: 2024 NORCE
# SPDX-License-Identifier: MIT

"""Test the parsing of .txt and .toml configuration files"""

import os
import sys
import filecmp
import pathlib
import subprocess

dirname: pathlib.Path = pathlib.Path(__file__).parent


def test_txt_toml():
    """See configs/spe11x_data_format.y (x in [a, b, c]; y in [txt, toml])"""
    if sys.version_info[1] < 11:
        print(
            "\nInput configuration files with toml extension requieres "
            + "a Python version of at least 3.11.\nTo run this test, "
            + "please update your Python version.\nYour Python version is "
            + f"3.{sys.version_info[1]}.{sys.version_info[2]}."
        )
        sys.exit()
    os.chdir(f"{dirname}/configs")
    for spe in ["spe11a", "spe11b", "spe11c"]:
        folder = []
        for ext in ["txt", "toml"]:
            subprocess.run(
                [
                    "pyopmspe11",
                    "-i",
                    f"{spe}_data_format.{ext}",
                    "-o",
                    f"{spe}_{ext}",
                    "-m",
                    "deck",
                ],
                check=True,
            )
            os.rename(
                f"{spe}_{ext}/deck/{spe.upper()}_{ext.upper()}.DATA",
                f"{spe}_{ext}/deck/{spe.upper()}.DATA",
            )
            folder.append(f"{dirname}/configs/{spe}_{ext}/deck")
        files = [f"{spe.upper()}.DATA", "TABLES.INC", "PERMX.INC"]
        if spe != "spe11a":
            files += ["PVBOUNDARIES.INC", "THCONR.INC"]
        if spe == "spe11c":
            files += ["GRID.INC"]
        match, mismatch, error = filecmp.cmpfiles(
            folder[0], folder[1], files, shallow=False
        )
        text = (
            f"Error for .txt and .tmol in {spe}: mismatch = {mismatch}; error = {error}"
        )
        assert match == files, text
