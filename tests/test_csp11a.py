# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the csp11a case"""

import os
from pyopmcsp11.core.pyopmcsp11 import main


def test_csp11a():
    """See configs/input.txt"""
    cwd = os.getcwd()
    os.chdir(f"{os.getcwd()}/tests/configs")
    main()
    os.chdir(cwd)
