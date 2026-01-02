# SPDX-FileCopyrightText: 2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Test the different models (immiscible, isothermal, convective, and complete)"""

import os
import pathlib
from mako.template import Template

testpth: pathlib.Path = pathlib.Path(__file__).parent


def test_models():
    """See configs/spe11x.mako"""
    if not os.path.exists(f"{testpth}/output"):
        os.mkdir(f"{testpth}/output")
    cases = ["a", "b", "c"]
    times = ["h", "y", "y"]
    models = ["immiscible", "isothermal", "convective", "complete"]
    flags = "-f 0 -g all -m deck_flow_data -t 1 -r 20,1,12"
    for case, time in zip(cases, times):
        for model in models:
            mytemplate = Template(filename=f"{testpth}/configs/spe11{case}.mako")
            var = {"model": model}
            filledtemplate = mytemplate.render(**var)
            name = f"{testpth}/output/spe11{case}_{model}"
            with open(
                f"{name}.toml",
                "w",
                encoding="utf8",
            ) as file:
                file.write(filledtemplate)
            os.system(
                f"pyopmspe11 -i {name}.toml -o {testpth}/output/spe11{case}_{model} {flags}"
            )
            data = f"{name}/spe11{case}_performance_spatial_map_1{time}.csv"
            assert os.path.exists(data)
