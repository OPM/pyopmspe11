# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the csp11c case"""

import os
import subprocess


def test_spe11b():
    """See configs/csp11c.txt"""
    cwd = os.getcwd()
    os.chdir(f"{os.getcwd()}/tests/configs")
    subprocess.run(["pyopmcsp11", "-i", "csp11c.txt", "-o", "csp11c"], check=True)
    os.chdir(cwd)
