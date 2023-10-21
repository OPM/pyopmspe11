# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the spe11a case"""

import os
import subprocess


def test_spe11a():
    """See configs/spe11a.txt"""
    cwd = os.getcwd()
    os.chdir(f"{cwd}/tests/configs")
    subprocess.run(["pyopmspe11", "-i", "spe11a.txt", "-o", "spe11a"], check=True)
    assert os.path.exists(f"{cwd}/tests/configs/spe11a/flow/SPE11A.UNRST")
    os.chdir(cwd)
