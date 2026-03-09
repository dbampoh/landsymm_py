"""Soil interpolation and output (stub).

Pseudocode:
1) For each soil input file:
   a) Read to geoarray (target gridlist).
   b) Map missing cells (diagnostic).
   c) Inpaint NaNs per variable.
   d) Validate no NaNs/negatives.
   e) Save output table with soil vars.
"""

# TODO:
# - implement soil interpolation per MATLAB remap_interp_soil
# - match output precision and formatting
from __future__ import annotations

from typing import Sequence

import numpy as np

from landsymm.common.interpolation import inpaint_nans
from landsymm.common.lpjg_io import read2geoarray, save_table
from landsymm.common.mapping_tools import xz_to_yxz


def interpolate_soil(
    files_soil: Sequence[str],
    gridlist,
    out_dir: str,
    out_dir_figs: str,
    inpaint_method: int,
    remap_ver: str,
    out_width: int,
    delimiter: str,
    overwrite: bool,
    fancy: bool,
) -> None:
    """Interpolate soil datasets to gridlist (TODO)."""
    n_cells = len(gridlist["list2map"])

    for file_soil in files_soil:
        soil_filename = file_soil.split("/")[-1]
        base, ext = soil_filename.rsplit(".", 1)
        ext = "." + ext
        print(f"{base}{ext}:")

        soil = read2geoarray(
            file_soil,
            verboseIfNoMat=False,
            force_mat_save=False,
            force_mat_nosave=True,
            target=gridlist,
            fill_missing_cells_with_nan=True,
        )

        # diagnostic map omitted for now (see remap_map_missing)
        out_soil = {
            "varNames": soil["varNames"],
            "lonlats": gridlist["lonlats"],
            "list2map": gridlist["list2map"],
            "garr_xv": np.full((n_cells, len(soil["varNames"])), np.nan),
        }
        for v, name in enumerate(soil["varNames"]):
            print(f"    Interpolating soil {name}...")
            tmp0_yx = xz_to_yxz(soil["garr_xv"][:, v:v+1], gridlist["mask_YX"].shape, soil["list2map"])[:, :, 0]
            tmp1_yx = inpaint_nans(tmp0_yx, inpaint_method)
            out_soil["garr_xv"][:, v] = tmp1_yx.ravel(order="F")[np.asarray(gridlist["list2map"]) - 1]

        if np.isnan(out_soil["garr_xv"]).any():
            raise RuntimeError("NaN in out_soil.garr_xv")
        if (out_soil["garr_xv"] < 0).any():
            raise RuntimeError("Negative in out_soil.garr_xv")

        out_header = ["lon", "lat"] + list(out_soil["varNames"])
        out_file_soil = f"{out_dir}/{base}.remapv{remap_ver}{ext}"

        print("    Saving...")
        save_table(
            out_header,
            out_soil,
            out_file_soil,
            outPrec=3,
            outWidth=out_width,
            delimiter=delimiter,
            overwrite=overwrite,
            fancy=fancy,
            progress_step_pct=20,
            verbose=False,
            gzip_output=False,
        )

    print("Done with soil.")
