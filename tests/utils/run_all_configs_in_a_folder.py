# SPDX-FileCopyrightText: 2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Generate all OPM input decks submitted to the benchmark (run this inside the benchmark folder)"""

import os
import subprocess


def main():
    """Entry point"""
    for i in range(5):
        for x in ("a", "b", "c"):
            if x in ("b", "c") and i == 4:
                continue
            fol = f"spe11{x}/{i}"
            if os.path.isdir(fol):
                subprocess.run(["rm", "-rf", fol], check=True)
    folder = os.path.abspath(".")
    for root, _, files in os.walk(folder):
        for i, file in enumerate(files):
            if file.endswith(".py"):
                continue
            confi = os.path.join(root, file)
            subprocess.run(
                [
                    "pyopmspe11",
                    "-i",
                    confi,
                    "-o",
                    root + f"/{i}",
                    "-f",
                    "0",
                    "-m",
                    "deck",
                ],
                check=True,
            )


if __name__ == "__main__":
    main()
