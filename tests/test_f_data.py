# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the script to write the data as required in the benchmark"""

from pyopmspe11.visualization.data import main


def test_f_data(test_d_f_g_workdir, monkeypatch):
    """Generate benchmark data output"""
    out = test_d_f_g_workdir / "output"
    flow = out / "flow"
    data = out / "data"
    assert (flow / "OUTPUT.UNRST").exists(), "Please run test_d_spe11b first"
    data.mkdir(exist_ok=True)
    monkeypatch.chdir(test_d_f_g_workdir)
    main()
    assert (data / "spe11b_time_series.csv").exists()
