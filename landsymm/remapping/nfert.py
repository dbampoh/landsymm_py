"""N fertilization processing (stub).

Pseudocode:
1) Map LandSyMM crop names to AgGRID types.
2) For combined crops, compute weighted Nfert (combine_nfert).
3) For normal crops, read NetCDF Napprate.
4) Enforce RF/IR equality.
5) Convert kg/ha to kg/m2.
"""

# TODO:
# - implement AgGRID mapping and NetCDF reads
# - ensure kg/ha -> kg/m2 conversion parity
from __future__ import annotations

from typing import Sequence

import numpy as np
from netCDF4 import Dataset


def aggrid_type_to_landsymm(crop: str) -> str:
    """Map LandSyMM crop to AgGRID type (TODO)."""
    if crop in {"CerealsC3", "StarchyRoots", "FruitAndVeg", "Sugarbeet", "OilOther", "ExtraCrop"}:
        return "wheat"
    if crop in {"CerealsC4", "Miscanthus", "Sugarcane"}:
        return "maize"
    if crop in {"Oilcrops", "Pulses", "OilNfix"}:
        return "soybean"
    if crop == "Rice":
        return "rice"
    if crop == "Sugar":
        return "combined_sugars"
    if crop == "OilPalm":
        # Use coffee because, like coffee, oil palm is a perennial tree.
        return "coffee"
    raise RuntimeError(f"What crop from AgGRID Nfert inputs should I use for {crop}?")


def combine_nfert(
    combined_crop: str,
    croparea_in,
    nfert_dir: str,
    gridlist,
) -> np.ndarray:
    """Combine Nfert for combined crop types (TODO)."""
    if combined_crop == "combined_sugars":
        croparea_names = ["Sugarbeet", "Sugarcane"]
        nfert_names = ["sugarbeet", "sugarcane"]
    else:
        raise RuntimeError(f"Not recognized: {combined_crop}")

    n_cells = croparea_in["garr_xv"].shape[0]
    area_xv = np.full((n_cells, len(croparea_names)), np.nan)
    for c, name in enumerate(croparea_names):
        is_this_crop = [name in v for v in croparea_in["varNames"]]
        if sum(is_this_crop) != 2:
            raise RuntimeError(
                f"Expected 2 croparea_in.varNames containing {name}; found {sum(is_this_crop)}"
            )
        area_xv[:, c] = np.sum(croparea_in["garr_xv"][:, np.array(is_this_crop)], axis=1)

    area_total_x = np.sum(area_xv, axis=1)
    if np.isnan(area_total_x).any():
        raise RuntimeError("Unexpected NaN in crop area")

    area_total_x[area_total_x == 0] = 1
    weights_xv = area_xv / np.tile(area_total_x.reshape(-1, 1), (1, len(croparea_names)))

    in_x = np.zeros((n_cells,))
    for c, name in enumerate(nfert_names):
        file_n = f"{nfert_dir}/agmip_{name}_apprate_fill_NPK_0.5.nc4"
        with Dataset(file_n, "r") as ds:
            var = ds.variables["Napprate"][:]
            if np.ma.isMaskedArray(var):
                var = var.filled(np.nan)
            map_yx = np.flipud(np.asarray(var))
        in_x = in_x + weights_xv[:, c] * map_yx.ravel(order="F")[np.asarray(gridlist["list2map"]) - 1]
    return in_x
