"""GeoArray data model (stub)."""
#
# TODO:
# - add validation helpers (shape, year_list length)
# - add constructors for maps_yxv vs maps_yxvy
# - add convenience methods (copy, subset years, reorder vars)
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence
import numpy as np


@dataclass
class GeoArray:
    var_names: list[str]
    year_list: Sequence[int] | None
    maps_yxv: np.ndarray | None
    maps_yxvy: np.ndarray | None
    lonlats: np.ndarray | None
    list2map: np.ndarray | None
    mask_yx: np.ndarray | None
