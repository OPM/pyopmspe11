# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the spe11b case"""

import os
from pyopmspe11.core.pyopmspe11 import main


def test_spe11b():
    """See configs/input.txt"""
    cwd = os.getcwd()
    os.chdir(f"{os.getcwd()}/tests/configs")
    main()
    os.chdir(cwd)
    assert os.path.exists(f"{cwd}/tests/configs/output/flow/OUTPUT.UNRST")
