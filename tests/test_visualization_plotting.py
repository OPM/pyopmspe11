# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

"""Test the visualisation scrip"""

import os
from pyopmcsp11.visualization.plotting import main


def test_visualization_plotting():
    """See visualization/plotting.py"""
    cwd = os.getcwd()
    os.chdir(f"{os.getcwd()}/tests/configs")
    main()
    os.chdir(cwd)
