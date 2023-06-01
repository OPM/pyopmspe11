# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the csp11b case"""

import os
import subprocess


def test_csp11b():
    """See configs/csp11b.txt"""
    cwd = os.getcwd()
    os.chdir(f"{os.getcwd()}/tests/configs")
    subprocess.run(["pyopmcsp11", "-i", "csp11b.txt", "-o", "csp11b"], check=True)
    os.chdir(cwd)
