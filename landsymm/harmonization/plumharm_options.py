"""Configuration model for PLUM harmonization (stub)."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence


@dataclass
class PlumHarmConfig:
    this_dir: str
    geodata_dir: str
    plum_dirs: Sequence[str]
    harm_dirs: Sequence[str] | None
    remap_lu_file: str
    remap_cropf_file: str
    remap_nfert_file: str
    base_year: int
    year1: int
    yearN: int
    allow_unveg: bool = False
    conserv_tol_pct: float = 0.2
    conserv_tol_area: float = 1e3
    norm2extra: float = 0.177
    use_latest_plum_mgmt: bool = True
    inpaint_method: int = 4
    fix_tiny_negs_tol_m2: float = 1.0
    save_halfdeg_mat: bool = True
    save_2deg_mat: bool = True
    save_halfdeg_txt: bool = False
    save_2deg_txt: bool = False
    out_prec: int = 6
    out_width: int = 1
    delimiter: str = " "
    overwrite: bool = True
    fancy: bool = False
    debug_areas: bool = False
    debug_nfert: bool = False
    debug_irrig: bool = False
    debugIJ_2deg: Sequence[int] | None = None
    dbCrop: str = ""
    verbose: bool = True
    combine_crops: bool = False
    fruitveg_sugar_2oil: bool = False
