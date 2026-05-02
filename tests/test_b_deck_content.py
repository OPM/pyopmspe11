# SPDX-FileCopyrightText: 2024-2026 NORCE Research AS
# SPDX-License-Identifier: MIT
# pylint: disable=R0912, R0914

"""Compare Cartesian and corner-point files used in the benchmark"""

import pathlib
import subprocess
import filecmp
import shutil

testpth = pathlib.Path(__file__).parent
mainpth = pathlib.Path(__file__).parents[1]


def test_b_deck_content(tmp_path, monkeypatch):
    """Compare generated decks against references"""
    monkeypatch.chdir(tmp_path)
    cases = {
        "a": ["r2_Cart_1cm_capmax2500Pa", "r3_cp_1cmish_capmax2500Pa"],
        "b": ["r1_Cart_10m", "r2_cp_10mish"],
        "c": ["r1_Cart_50m-50m-10m", "r2_cp_50m-50m-8mish"],
    }
    for c, configs in cases.items():
        casef = tmp_path / f"spe11{c}"
        casef.mkdir()
        for case in configs:
            run = casef / case
            subprocess.run(
                [
                    "pyopmspe11",
                    "-i",
                    str(mainpth / "benchmark" / f"spe11{c}" / f"{case}.toml"),
                    "-o",
                    str(run),
                    "-m",
                    "deck",
                    "-f",
                    "0",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            ref = testpth / "decks" / f"spe11{c}" / case
            grid_case = ref / "GRID.INC"
            if case == "r2_cp_50m-50m-8mish" and not grid_case.is_file():
                grid_case_zip = ref / "GRID.INC.zip"
                shutil.unpack_archive(grid_case_zip, ref)
            files = sorted(
                p.name
                for p in ref.iterdir()
                if p.is_file() and not p.name.endswith(".zip")
            )
            if "PVBOUNDARIES.INC" in files:
                # Different Python versions/operative systems might results in
                # floating differences, then compare the sum up to a tolerance
                # and the indices
                files.remove("PVBOUNDARIES.INC")

                ref_porv, ref_ind = [], []
                with open(ref / "PVBOUNDARIES.INC", "r", encoding="utf8") as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("--"):
                            continue
                        if line.startswith("PORV"):
                            parts = line.split()
                            ref_porv.append(float(parts[1]))
                            ref_ind.append(" ".join(parts[2:]))

                test_porv, test_ind = [], []
                with open(run / "PVBOUNDARIES.INC", "r", encoding="utf8") as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("--"):
                            continue
                        if line.startswith("PORV"):
                            parts = line.split()
                            test_porv.append(float(parts[1]))
                            test_ind.append(" ".join(parts[2:]))

                assert len(ref_porv) == len(test_porv)
                for rpv, tpv, rind, tind in zip(ref_porv, test_porv, ref_ind, test_ind):
                    assert abs(rpv - tpv) < 1e-6
                    assert rind == tind

            _, mismatch, error = filecmp.cmpfiles(run, ref, files, shallow=False)

            assert not mismatch and not error, (
                f"Config mismatch for spe11{c}/{case}\n"
                f"Mismatched contents: {mismatch}\n"
                f"Errors: {error}"
            )
