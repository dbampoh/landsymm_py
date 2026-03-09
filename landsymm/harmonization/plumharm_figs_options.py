"""Configuration model for PLUMharmFigs."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence


@dataclass
class PlumHarmFigsConfig:
    # Needed from PLUMharm_options.m
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
    harms_figs_dir: str
    combine_crops: bool = False
    fruitveg_sugar_2oil: bool = False
    allow_unveg: bool = False
    norm2extra: float = 0.177
    timeseries_legend_loc: str = "best"
    runlist_legend: Sequence[str] | None = None
    save_geotiffs: bool = False
    this_ver: str = ""
    three_years: Sequence[int] | None = None
    year_list_baseline_to_plot: Sequence[int] | None = None
