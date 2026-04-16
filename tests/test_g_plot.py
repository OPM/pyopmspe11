# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the script to plot the data"""

from pyopmspe11.visualization.plotting import main


def test_g_plot(test_d_f_g_workdir, monkeypatch):
    """Generate benchmark plots"""
    out = test_d_f_g_workdir / "output"
    data = out / "data"
    figs = out / "figures"
    assert data.exists(), "Please run test_f_data first"
    figs.mkdir(exist_ok=True)
    monkeypatch.chdir(test_d_f_g_workdir)
    main()
    assert (figs / "spe11b_sparse_data.png").exists()
