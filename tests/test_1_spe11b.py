# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the spe11b case"""

import os
import pathlib
from pyopmspe11.core.pyopmspe11 import main

dirname: pathlib.Path = pathlib.Path(__file__).parent


def test_spe11b():
    """See configs/input.txt"""
    os.chdir(f"{dirname}/configs")
    main()
    assert os.path.exists(f"{dirname}/configs/output/flow/OUTPUT.UNRST")
