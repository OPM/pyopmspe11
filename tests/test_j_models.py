# SPDX-FileCopyrightText: 2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the different models (immiscible, isothermal, convective, and complete)"""

import pathlib
from mako.template import Template

from pyopmspe11.core.pyopmspe11 import main

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_j_models(tmp_path, monkeypatch):
    """Run different physical models and validate outputs"""
    cases = ("a", "b", "c")
    times = ("h", "y", "y")
    models = ("immiscible", "isothermal", "convective", "complete")
    flags = ["-f", "0", "-g", "all", "-m", "deck_flow_data", "-t", "1", "-r", "20,1,12"]
    for case, time in zip(cases, times):
        template = Template(filename=f"{testpth}/configs/spe11{case}.mako")
        for model in models:
            monkeypatch.chdir(tmp_path)
            cfg = tmp_path / f"spe11{case}_{model}.toml"
            cfg.write_text(template.render(model=model), encoding="utf8")
            out = tmp_path / f"spe11{case}_{model}"
            main(["-i", str(cfg), "-o", str(out.name)] + flags)
            data = out / f"spe11{case}_performance_spatial_map_1{time}.csv"
            assert data.exists()
