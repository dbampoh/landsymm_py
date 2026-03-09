"""Gridlist utilities (stub)."""
#
# TODO:
# - parse gridlist text file into lonlats/list2map/mask
# - support lat_extent and orientation flags
# - ensure list2map ordering matches MATLAB
from __future__ import annotations


def load_gridlist(
    path: str,
    *,
    xres: float | None = None,
    yres: float | None = None,
    lat_orient: str = "",
    lon_orient: str = "",
    drop_northpole: bool = False,
    drop_southpole: bool = False,
    lons_centered_on_180: bool = False,
    verbose: bool = False,
):
    """Load gridlist and list2map/mask with MATLAB-like indexing."""
    from .lpjg_io import read2geoarray

    out = read2geoarray(
        path,
        xres=xres if xres is not None else float("nan"),
        yres=yres if yres is not None else float("nan"),
        lat_orient=lat_orient,
        lon_orient=lon_orient,
        verbose=verbose,
        force_mat_save=True,
        force_mat_nosave=True,
        drop_northpole=drop_northpole,
        drop_southpole=drop_southpole,
        lons_centered_on_180=lons_centered_on_180,
        gridlist_warn=False,
    )
    return out
