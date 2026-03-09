"""Configuration model for PLUMharm2LPJG."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence


@dataclass
class PlumHarm2LPJGConfig:
    # Needed from PLUMharm_options.m
    this_dir: str
    plum_dirs: Sequence[str]
    harm_dirs: Sequence[str] | None
    base_year: int
    year1: int
    yearN: int
    combine_crops: bool = False
    fruitveg_sugar_2oil: bool = False
    out_prec: int = 6
    out_width: int = 1
    delimiter: str = " "
    overwrite: bool = True
    fancy: bool = False

    # PLUMharm2LPJG-specific options
    do_gzip: bool = False
    donation_order: Sequence[str] = ("PASTURE", "NATURAL", "BARREN")
    save_every_pct: int = 1
    someofall: bool = True
    forLPJG_dirs: Sequence[str] | None = None
    verbose_write: bool = False
    yStep: int = 1
    y1_pre: int | None = None
