"""Configuration model for remapping (stub)."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence


@dataclass
class RemapConfig:
    geodata_dir: str
    year_list_out: Sequence[int]
    lu_source: str
    out_dir_top: str
    fill_unveg: float
    interp_test_y1: int | None = None
    interp_test_yN: int | None = None
    inpaint_method: int = 4
    force_all_rainfed: bool = False
    plum_setaside_frac: float = 0.0
    remap_ver: str = ""
    this_ver: str = ""
    file_gridlist_out: str = ""
    file_gridlist_climate: str | None = None
    files_soil: Sequence[str] = ()
    out_prec: int = 6
    out_width: int = 1
    delimiter: str = " "
    overwrite: bool = True
    fancy: bool = False
    do_gzip: bool = False
    max_nan_frac: float = 0.01
    log_file: str | None = None
