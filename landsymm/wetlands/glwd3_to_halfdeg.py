"""Aggregate GLWD3 30-arcsec raster to half-degree wetland fractions.

Why this module exists (OPTIONAL stage)
=======================================
This is the first step of **Stage 4 of the landsymm_py pipeline, which
is optional**. Run it only if your downstream LPJ-GUESS runs need an
explicit PEATLAND land-cover class — typically because you are running
coupled LPJ-GUESS ↔ IMOGEN climate simulations where peatland CH₄
emissions feed back into the IMOGEN climate model.

Source: Lehner & Döll (2004), J. Hydrology 296:1-22 (GLWD3).

Python port of glwd3_to_halfdeg.R.

Reads the GLWD3 Level-3 GeoTIFF (12 land-cover classes at ~1 km resolution),
aggregates each class to 0.5-degree cell fractions, then applies the Lehner (2004)
scale factors to produce total wetland fraction per half-degree cell.

Two wetland products are generated (mirroring the R script):
  - "with forests": includes Swamp Forest (class 5) and Coastal Wetland (class 6)
  - "without forests": excludes those forest-type wetlands

The "with forests" product is appropriate when peatland is being carved from
NATURAL as a distinct land type (forest wetlands are genuine wetlands whose tree
cover is handled by the vegetation dynamics model). The "without forests" product
avoids potential double-counting if the model already explicitly represents swamp
forests and mangroves within its natural vegetation PFTs.

Output is saved as a single NetCDF file (peatland_halfdeg.nc) containing:
  - 12 layers of raw per-class fractions
  - 'wetland_frac' (with forests) and 'wetland_frac_noforests' (without forests)
"""
from __future__ import annotations

import argparse
import os
import time

import numpy as np

GLWD3_CLASSES = {
    1: "Lake",
    2: "Reservoir",
    3: "River",
    4: "Freshwater Marsh, Floodplain",
    5: "Swamp Forest, Flooded Forest",
    6: "Coastal Wetland",
    7: "Pan, Brackish/Saline Wetland",
    8: "Bog, Fen, Mire (Peatland)",
    9: "Intermittent Wetland/Lake",
    10: "50-100% Wetland",
    11: "25-50% Wetland",
    12: "Wetland Complex (0-25%)",
}

# Lehner (2004) Table 5 — wetland fraction INCLUDING forest wetland types
SCALE_FACTORS_WFORESTS = np.array([
    0.0,    # 1  Lake
    0.0,    # 2  Reservoir
    0.0,    # 3  River
    1.0,    # 4  Freshwater Marsh, Floodplain
    1.0,    # 5  Swamp Forest, Flooded Forest
    1.0,    # 6  Coastal Wetland
    1.0,    # 7  Pan, Brackish/Saline Wetland
    1.0,    # 8  Bog, Fen, Mire (Peatland)
    1.0,    # 9  Intermittent Wetland/Lake
    0.75,   # 10 50-100% Wetland
    0.375,  # 11 25-50% Wetland
    0.125,  # 12 Wetland Complex (0-25%)
])

# Same scale factors but EXCLUDING forest wetland types (classes 5 & 6)
SCALE_FACTORS_NOFORESTS = np.array([
    0.0,    # 1  Lake
    0.0,    # 2  Reservoir
    0.0,    # 3  River
    1.0,    # 4  Freshwater Marsh, Floodplain
    0.0,    # 5  Swamp Forest, Flooded Forest  — EXCLUDED
    0.0,    # 6  Coastal Wetland               — EXCLUDED
    1.0,    # 7  Pan, Brackish/Saline Wetland
    1.0,    # 8  Bog, Fen, Mire (Peatland)
    1.0,    # 9  Intermittent Wetland/Lake
    0.75,   # 10 50-100% Wetland
    0.375,  # 11 25-50% Wetland
    0.125,  # 12 Wetland Complex (0-25%)
])


def aggregate_glwd3_to_halfdeg(
    glwd3_path: str,
    output_dir: str,
    verbose: bool = True,
) -> str:
    """Aggregate GLWD3 GeoTIFF to half-degree wetland fractions.

    Parameters
    ----------
    glwd3_path : str
        Path to the GLWD3 GeoTIFF file (glwd_3.tif).
    output_dir : str
        Directory for output files.
    verbose : bool
        Print progress messages.

    Returns
    -------
    str
        Path to the output NetCDF file.
    """
    import rasterio
    from netCDF4 import Dataset

    t0 = time.perf_counter()

    if verbose:
        print(f"Reading GLWD3 raster: {glwd3_path}")
    with rasterio.open(glwd3_path) as src:
        data = src.read(1)
        nodata = src.nodata
    nrows_native, ncols_native = data.shape
    if verbose:
        print(f"  Native resolution: {nrows_native} x {ncols_native}")

    # The GLWD3 GeoTIFF may be 21600x43200 (30-arcsec global)
    # or it could have slight variations. Compute the aggregation factor.
    target_nrows = 360  # half-degree: 180 / 0.5
    target_ncols = 720  # half-degree: 360 / 0.5
    ag_row = nrows_native // target_nrows
    ag_col = ncols_native // target_ncols

    if ag_row * target_nrows != nrows_native or ag_col * target_ncols != ncols_native:
        raise RuntimeError(
            f"GLWD3 dimensions ({nrows_native}x{ncols_native}) are not evenly "
            f"divisible by half-degree target ({target_nrows}x{target_ncols})"
        )

    n_per_cell = ag_row * ag_col
    if verbose:
        print(f"  Aggregation factor: {ag_row} x {ag_col} = {n_per_cell} pixels per half-deg cell")

    # Replace nodata with 0 (not a valid class)
    if nodata is not None:
        data[data == int(nodata)] = 0

    # Aggregate each class to half-degree fraction
    halfdeg_fracs = np.zeros((12, target_nrows, target_ncols), dtype=np.float64)

    for cls_idx in range(12):
        cls_id = cls_idx + 1
        if verbose:
            print(f"  Aggregating class {cls_id}: {GLWD3_CLASSES[cls_id]}")
        binary = (data == cls_id).astype(np.float64)
        reshaped = binary.reshape(target_nrows, ag_row, target_ncols, ag_col)
        halfdeg_fracs[cls_idx] = reshaped.sum(axis=(1, 3)) / n_per_cell

    # Apply scale factors to get wetland fractions — both products
    scaled_wf = halfdeg_fracs * SCALE_FACTORS_WFORESTS[:, None, None]
    wetland_frac = scaled_wf.sum(axis=0)

    scaled_nf = halfdeg_fracs * SCALE_FACTORS_NOFORESTS[:, None, None]
    wetland_frac_noforests = scaled_nf.sum(axis=0)

    if verbose:
        print(f"  With forests — wetland fraction range: [{wetland_frac.min():.6f}, {wetland_frac.max():.6f}]")
        print(f"    Non-zero cells: {np.sum(wetland_frac > 0):,}")
        print(f"  Without forests — wetland fraction range: [{wetland_frac_noforests.min():.6f}, {wetland_frac_noforests.max():.6f}]")
        print(f"    Non-zero cells: {np.sum(wetland_frac_noforests > 0):,}")

    # Save as NetCDF
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "peatland_halfdeg.nc")

    if verbose:
        print(f"  Saving to {out_path}")

    lons = np.arange(-179.75, 180.0, 0.5)
    lats = np.arange(89.75, -90.0, -0.5)

    with Dataset(out_path, "w", format="NETCDF4") as nc:
        nc.createDimension("lat", target_nrows)
        nc.createDimension("lon", target_ncols)
        nc.createDimension("class", 12)

        lat_var = nc.createVariable("lat", "f8", ("lat",))
        lat_var[:] = lats
        lat_var.units = "degrees_north"
        lat_var.long_name = "latitude"

        lon_var = nc.createVariable("lon", "f8", ("lon",))
        lon_var[:] = lons
        lon_var.units = "degrees_east"
        lon_var.long_name = "longitude"

        cls_var = nc.createVariable("class_id", "i4", ("class",))
        cls_var[:] = np.arange(1, 13)
        cls_var.long_name = "GLWD3 class identifier"

        frac_var = nc.createVariable(
            "class_frac", "f8", ("class", "lat", "lon"),
            zlib=True, complevel=4,
        )
        frac_var[:] = halfdeg_fracs
        frac_var.long_name = "Raw fraction of half-degree cell in each GLWD3 class"
        frac_var.units = "fraction"

        scaled_wf_var = nc.createVariable(
            "scaled_frac_wforests", "f8", ("class", "lat", "lon"),
            zlib=True, complevel=4,
        )
        scaled_wf_var[:] = scaled_wf
        scaled_wf_var.long_name = "Scaled wetland fraction per class (with forest types)"
        scaled_wf_var.units = "fraction"

        scaled_nf_var = nc.createVariable(
            "scaled_frac_noforests", "f8", ("class", "lat", "lon"),
            zlib=True, complevel=4,
        )
        scaled_nf_var[:] = scaled_nf
        scaled_nf_var.long_name = "Scaled wetland fraction per class (without forest types)"
        scaled_nf_var.units = "fraction"

        wf_var = nc.createVariable(
            "wetland_frac", "f8", ("lat", "lon"),
            zlib=True, complevel=4,
        )
        wf_var[:] = wetland_frac
        wf_var.long_name = "Total wetland fraction (with forest types: classes 4-12)"
        wf_var.units = "fraction"

        nf_var = nc.createVariable(
            "wetland_frac_noforests", "f8", ("lat", "lon"),
            zlib=True, complevel=4,
        )
        nf_var[:] = wetland_frac_noforests
        nf_var.long_name = "Total wetland fraction (without forest types: excludes classes 5,6)"
        nf_var.units = "fraction"

        nc.title = "GLWD3 wetland fractions aggregated to 0.5-degree resolution"
        nc.source = "GLWD Level 3 (Lehner & Döll 2004)"
        nc.history = "Aggregated from 30-arcsec to 0.5-deg with Lehner (2004) Table 5 scale factors"

        for i, (cls_id, name) in enumerate(GLWD3_CLASSES.items()):
            nc.setncattr(f"class_{cls_id}_name", name)
            nc.setncattr(f"class_{cls_id}_scale_wforests", float(SCALE_FACTORS_WFORESTS[i]))
            nc.setncattr(f"class_{cls_id}_scale_noforests", float(SCALE_FACTORS_NOFORESTS[i]))

    elapsed = time.perf_counter() - t0
    if verbose:
        print(f"  Done ({elapsed:.1f}s)")

    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Aggregate GLWD3 to half-degree wetland fractions."
    )
    parser.add_argument(
        "--glwd3-path",
        default=None,
        help="Path to GLWD3 GeoTIFF (default: data/geodata_py/glwd3/glwd_3.tif)",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: data/geodata_py/glwd3/)",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    from landsymm.config import get_geodata_dir

    glwd3_dir = str(get_geodata_dir() / "glwd3")
    glwd3_path = args.glwd3_path or str(get_geodata_dir() / "glwd3" / "glwd_3.tif")
    output_dir = args.output_dir or glwd3_dir

    aggregate_glwd3_to_halfdeg(glwd3_path, output_dir, verbose=not args.quiet)


if __name__ == "__main__":
    main()
