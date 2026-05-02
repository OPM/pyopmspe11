#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2023-2026 NORCE Research AS
# SPDX-License-Identifier: GPL-3.0
# pylint: disable=C0103, R0902

"""Central configuration models for pyopmspe11"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional


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
    dims: List[float]
    x_n: List[int]
    y_n: List[int]
    z_n: List[int]
    temperature: List[float]
    datum: float
    pressure: float
    kzMult: float
    diffusion: List[float]
    dispersion: List[float]
    radius: List[float]
    wellCoord: List[List[float]]
    krw: str
    krn: str
    pcap: str
    s_w: str
    safu: List[List[float]]
    rock: List[List[float]]
    inj: List[List[float]]
    # ------------------------------------------------------------------
    # TOML configuration optional (e.g., bc spe11a, convective mixing)
    # ------------------------------------------------------------------
    spe11aBC: Optional[float] = 0
    drsdtcon: Optional[List[List[str]]] = None
    elevation: Optional[float] = None
    backElevation: Optional[float] = None
    rockCond: Optional[List[float]] = None
    widthBuffer: Optional[float] = None
    rockExtra: Optional[List[float]] = None
    pvAdded: Optional[float] = None
    wellCoordF: Optional[List[List[float]]] = None
    # ------------------------------------------------------------------
    # SPE11 geometry and observation setup
    # ------------------------------------------------------------------
    maxelevation: float = 0
    cut: Optional[float] = 0
    nxyz: List[int] = field(default_factory=lambda: [0, 0, 0])
    boxa: List[List[float]] = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])
    boxb: List[List[float]] = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])
    boxc: List[List[float]] = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])
    sensors: List[List[float]] = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])
    sensorijk: List[List[int]] = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])
    wellijk: List[List[int]] = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])
    wellijkf: List[List[int]] = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])
    wellkh: Optional[List[int]] = field(default_factory=lambda: [])
    # ------------------------------------------------------------------
    # Miscellaneous runtime flags and metadata
    # ------------------------------------------------------------------
    pat: Path = Path(__file__).resolve().parents[1]  # Do not overwritte
    tuning: bool = False
    deckfol: str = "output"
    compact_dx: bool = False
