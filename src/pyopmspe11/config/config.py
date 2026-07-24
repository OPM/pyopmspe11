# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: GPL-3.0
# pylint: disable=C0103, R0902

"""Central configuration models for pyopmspe11"""

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(slots=True)
class Config:
    """Combines CLI options, TOML inputs, and derived runtime settings"""

    # ------------------------------------------------------------------
    # CLI configuration (normalized)
    # ------------------------------------------------------------------
    fol: str
    generate: str
    mode: str
    resolution: str
    time_data: str
    dt_data: float
    lower: bool
    subfolders: str
    # ------------------------------------------------------------------
    # TOML configuration (simulation setup)
    # ------------------------------------------------------------------
    flow: str
    spe11: str
    version: str
    model: str
    grid: str
    dims: list[float]
    x_n: list[int]
    y_n: list[int]
    z_n: list[int]
    temperature: list[float]
    datum: float
    pressure: float
    kzMult: float
    diffusion: list[float]
    dispersion: list[float]
    radius: list[float]
    wellCoord: list[list[float]]
    krw: str
    krn: str
    pcap: str
    s_w: str
    safu: list[list[float]]
    rock: list[list[float]]
    inj: list[list[float]]
    # ------------------------------------------------------------------
    # TOML configuration optional (e.g., bc spe11a, convective mixing)
    # ------------------------------------------------------------------
    spe11aBC: float | None = 0
    drsdtcon: list[list[str]] | None = None
    elevation: float | None = None
    backElevation: float | None = None
    rockCond: list[float] | None = None
    widthBuffer: float | None = None
    rockExtra: list[float] | None = None
    pvAdded: float | None = None
    wellCoordF: list[list[float]] | None = None
    # ------------------------------------------------------------------
    # SPE11 geometry and observation setup
    # ------------------------------------------------------------------
    maxelevation: float = 0
    cut: float | None = 0
    nxyz: list[int] = field(default_factory=lambda: [0, 0, 0])
    boxa: list[list[float]] = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])
    boxb: list[list[float]] = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])
    boxc: list[list[float]] = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])
    sensors: list[list[float]] = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])
    sensorijk: list[list[int]] = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])
    wellijk: list[list[int]] = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])
    wellijkf: list[list[int]] = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])
    wellkh: list[int] | None = field(default_factory=list)
    # ------------------------------------------------------------------
    # Miscellaneous runtime flags and metadata
    # ------------------------------------------------------------------
    pat: Path = Path(__file__).resolve().parents[1]  # Do not overwritte
    tuning: bool = False
    deckfol: str = "output"
    compact_dx: bool = False
