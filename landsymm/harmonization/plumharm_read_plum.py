"""PLUM input reading and caching (MATLAB parity)."""
from __future__ import annotations

from typing import Any, Dict, Sequence, Tuple

from .plumharm_pp_read_plum import pp_read_plum


def read_plum_inputs(
    in_dir: str,
    base_year: int,
    year_list,
    land_area_yx,
    lu_names,
    plum_to_lpjg,
    lpjg_crops,
    is_2deg: bool,
    bare_frac_y0_yx,
    norm2extra: float,
    inpaint_method: int,
    is_orig: bool,
    fruitveg_sugar_2oil: bool,
    allow_unveg: bool,
):
    """Read PLUM LandCoverFract/CropFract/Fert/Irrig with caching."""
    return pp_read_plum(
        in_dir,
        base_year,
        year_list,
        land_area_yx,
        lu_names,
        plum_to_lpjg,
        lpjg_crops,
        is_2deg,
        bare_frac_y0_yx,
        norm2extra,
        inpaint_method,
        "" if is_orig else "harm.",
        is_orig,
        fruitveg_sugar_2oil,
        allow_unveg,
    )
