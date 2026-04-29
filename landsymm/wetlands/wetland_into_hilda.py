"""Insert GLWD3 wetland/peatland fraction into HILDA+ remap land-use data.

Why this module exists (OPTIONAL stage)
=======================================
This module belongs to **Stage 4 of the landsymm_py pipeline, which is
optional**. Run it only if your downstream LPJ-GUESS runs need an
explicit PEATLAND land-cover class — typically because you are running
**coupled LPJ-GUESS ↔ IMOGEN climate simulations** in which peatland
CH₄ emissions feed back into the IMOGEN intermediate-complexity climate
model alongside CO₂ and N₂O, and that modified climate then drives
subsequent LPJ-GUESS ecosystem responses (Rabin et al., 2020,
Earth Syst. Dynam. 11:357-376; Wania et al., 2009, 2010 for the
peatland CH₄ submodel in LPJ-GUESS).

If your LPJ-GUESS runs use prescribed (offline) climate forcing and
do not need explicit peatland CH₄ accounting, you can skip this stage —
the Stage-2 HILDA+ remap output is already a complete, LPJ-GUESS-ready
historical land-use file.

What this module does
=====================
HILDA+ does not distinguish peatland from other natural vegetation, so
this module carves a PEATLAND fraction out of NATURAL using GLWD3
(Lehner & Döll, 2004) as the source of peatland extent. The
companion ``wetland_into_forLPJG.py`` performs the same operation on
the harmonized PLUM scenario landcover files (Stage 3 outputs).

Python port of glwd3_wetland_into_luh2.R — Approach H only,
adapted for HILDA+ remap LU format (NATURAL, CROPLAND, PASTURE, BARREN).

Approach H: Time-invariant peatland carved from NATURAL.
  1. For each gridcell (Lon, Lat), find the minimum NATURAL fraction
     across ALL years in the dataset.
  2. Set PEATLAND = min(min_NATURAL, GLWD3_wetland_frac) — guaranteed
     to always be available from NATURAL.
  3. NATURAL_new = NATURAL_old - PEATLAND at every year.
  4. CROPLAND, PASTURE, BARREN are unchanged.
  5. Fractions sum to 1.0 at every gridcell and year.

The min-over-years approach guarantees that PEATLAND is bounded by the
minimum natural-vegetation availability and never exceeds it in any
year, even under scenarios where natural vegetation declines over time.

Output is written alongside the original LU file with a distinguishing
``_peatland`` suffix to keep peatland and non-peatland baseline files
separate so a single project can support both kinds of LPJ-GUESS runs
side by side.
"""
from __future__ import annotations

import argparse
import os
import time

import numpy as np
import pandas as pd


WETLAND_VARS = {
    "wforests": "wetland_frac",
    "noforests": "wetland_frac_noforests",
}


def _load_wetland_frac_from_nc(
    nc_path: str,
    wetland_product: str = "wforests",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load the half-degree wetland fraction grid from the aggregated NetCDF.

    Parameters
    ----------
    nc_path : str
        Path to the NetCDF file produced by glwd3_to_halfdeg.py.
    wetland_product : str
        Which wetland product to use: "wforests" (includes forest wetland
        types) or "noforests" (excludes classes 5 & 6).

    Returns (wetland_frac_2d, lats_1d, lons_1d).
    wetland_frac_2d is (360, 720) with lats descending from 89.75 to -89.75.
    """
    from netCDF4 import Dataset

    var_name = WETLAND_VARS.get(wetland_product)
    if var_name is None:
        raise ValueError(
            f"Unknown wetland_product '{wetland_product}'. "
            f"Choose from: {list(WETLAND_VARS.keys())}"
        )

    with Dataset(nc_path, "r") as nc:
        if var_name not in nc.variables:
            # Backward compatibility: older files may only have "wetland_frac"
            if wetland_product == "wforests" and "wetland_frac" in nc.variables:
                var_name = "wetland_frac"
            else:
                raise RuntimeError(
                    f"Variable '{var_name}' not found in {nc_path}. "
                    f"Available: {list(nc.variables.keys())}"
                )
        wf = np.asarray(nc.variables[var_name][:], dtype=np.float64)
        lats = np.asarray(nc.variables["lat"][:], dtype=np.float64)
        lons = np.asarray(nc.variables["lon"][:], dtype=np.float64)
    return wf, lats, lons


def _extract_wetland_frac_for_coords(
    wf_2d: np.ndarray,
    lats_grid: np.ndarray,
    lons_grid: np.ndarray,
    lon_vals: np.ndarray,
    lat_vals: np.ndarray,
) -> np.ndarray:
    """Extract wetland fraction for each (Lon, Lat) coordinate pair.

    The NetCDF grid has cell centers at -179.75, -179.25, ..., 179.75 (lon)
    and 89.75, 89.25, ..., -89.75 (lat, descending).
    """
    lon_idx = np.searchsorted(lons_grid, lon_vals)
    lon_idx = np.clip(lon_idx, 0, len(lons_grid) - 1)

    lat_idx = np.searchsorted(-lats_grid, -lat_vals)
    lat_idx = np.clip(lat_idx, 0, len(lats_grid) - 1)

    return wf_2d[lat_idx, lon_idx]


def insert_wetland_approach_h(
    lu_path: str,
    wetland_nc_path: str,
    output_path: str,
    wetland_product: str = "wforests",
    verbose: bool = True,
) -> str:
    """Insert GLWD3 wetland into HILDA+ remap LU using Approach H.

    Parameters
    ----------
    lu_path : str
        Path to the original HILDA+ remap LU file (space-delimited,
        columns: Lon Lat Year NATURAL CROPLAND PASTURE BARREN).
    wetland_nc_path : str
        Path to the peatland_halfdeg.nc file from Step 1.
    output_path : str
        Path for the output LU file with PEATLAND column added.
    wetland_product : str
        Which wetland product to use: "wforests" (includes forest wetland
        types) or "noforests" (excludes classes 5 & 6).
    verbose : bool
        Print progress.

    Returns
    -------
    str
        Path to the output file.
    """
    t0 = time.perf_counter()

    if verbose:
        print(f"Loading wetland fraction grid: {wetland_nc_path}")
        print(f"  Wetland product: {wetland_product}")
    wf_2d, lats_grid, lons_grid = _load_wetland_frac_from_nc(
        wetland_nc_path, wetland_product=wetland_product
    )

    if verbose:
        print(f"Reading HILDA+ remap LU: {lu_path}")
    lu = pd.read_csv(lu_path, sep=r"\s+")

    expected_cols = {"Lon", "Lat", "Year", "NATURAL", "CROPLAND", "PASTURE", "BARREN"}
    if not expected_cols.issubset(set(lu.columns)):
        missing = expected_cols - set(lu.columns)
        raise RuntimeError(f"Missing columns in LU file: {missing}")

    n_rows = len(lu)
    years = sorted(lu["Year"].unique())
    n_years = len(years)
    if verbose:
        print(f"  {n_rows:,} rows, {n_years} years ({years[0]}-{years[-1]})")

    # Extract the GLWD3 wetland fraction for each coordinate pair
    # (same fraction for all years since GLWD3 is time-invariant)
    if verbose:
        print("Extracting GLWD3 wetland fraction per gridcell...")
    peatland_frac_raw = _extract_wetland_frac_for_coords(
        wf_2d, lats_grid, lons_grid,
        lu["Lon"].values, lu["Lat"].values,
    )

    # --- Approach H ---
    # Step 1: For each gridcell, find minimum NATURAL across all years
    if verbose:
        print("Computing minimum NATURAL per gridcell across all years...")
    min_nat = lu.groupby(["Lon", "Lat"])["NATURAL"].transform("min")

    # Step 2: PEATLAND = min(min_NATURAL, GLWD3_frac)
    peatland = np.minimum(min_nat.values, peatland_frac_raw)

    # Sanity: ensure non-negative
    peatland = np.maximum(peatland, 0.0)

    # Step 3: NATURAL_new = NATURAL - PEATLAND
    new_natural = lu["NATURAL"].values - peatland

    # Check for negatives (should not happen by construction, but verify)
    n_neg = np.sum(new_natural < 0)
    if n_neg > 0:
        if verbose:
            print(f"  WARNING: {n_neg} negative NATURAL values detected, clamping to 0")
        neg_mask = new_natural < 0
        peatland[neg_mask] = lu["NATURAL"].values[neg_mask]
        new_natural[neg_mask] = 0.0

    # Step 4: Build output dataframe
    lu_out = lu.copy()
    lu_out["NATURAL"] = new_natural
    lu_out["PEATLAND"] = peatland

    # Reorder columns: Lon Lat Year NATURAL CROPLAND PASTURE PEATLAND BARREN
    col_order = ["Lon", "Lat", "Year", "NATURAL", "CROPLAND", "PASTURE", "PEATLAND", "BARREN"]
    lu_out = lu_out[col_order]

    # Verify fractions sum to 1.0
    frac_cols = ["NATURAL", "CROPLAND", "PASTURE", "PEATLAND", "BARREN"]
    row_sums = lu_out[frac_cols].sum(axis=1)
    max_dev = np.max(np.abs(row_sums - 1.0))
    if verbose:
        print(f"  Max deviation from sum=1.0: {max_dev:.2e}")
    if max_dev > 1e-10:
        print(f"  WARNING: Row sums deviate from 1.0 by up to {max_dev:.6e}")

    # Summary statistics
    if verbose:
        n_with_peatland = np.sum(peatland > 0)
        n_cells_with_peatland = len(lu_out.loc[lu_out["PEATLAND"] > 0, ["Lon", "Lat"]].drop_duplicates())
        total_rows = len(lu_out)
        print(f"  Rows with PEATLAND > 0: {n_with_peatland:,} / {total_rows:,}")
        print(f"  Unique gridcells with PEATLAND > 0: {n_cells_with_peatland:,}")
        print(f"  PEATLAND range: [{peatland[peatland > 0].min():.6f}, {peatland.max():.6f}]" if n_with_peatland > 0 else "  No peatland assigned")

    # Write output
    if verbose:
        print(f"Writing output: {output_path}")

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    with open(output_path, "w") as f:
        header = " ".join(f"{c:>10s}" for c in col_order)
        f.write(header + "\n")
        for i in range(0, n_rows, 50000):
            chunk = lu_out.iloc[i:i + 50000]
            lines = []
            for _, row in chunk.iterrows():
                parts = [
                    f"{row['Lon']:10.2f}",
                    f"{row['Lat']:10.2f}",
                    f"{int(row['Year']):10d}",
                ]
                for col in frac_cols:
                    parts.append(f"{row[col]:10.6f}")
                lines.append(" ".join(parts))
            f.write("\n".join(lines) + "\n")
            if verbose and i > 0 and i % 500000 == 0:
                print(f"  Written {i:,} / {n_rows:,} rows...")

    elapsed = time.perf_counter() - t0
    if verbose:
        print(f"Done ({elapsed:.1f}s)")

    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Insert GLWD3 wetland/peatland into HILDA+ remap LU data (Approach H). "
            "Output file is placed alongside the original with a '_peatland' tag."
        )
    )
    parser.add_argument(
        "--lu-path",
        default=None,
        help="Path to HILDA+ remap LU file (default: auto-detect in data/)",
    )
    parser.add_argument(
        "--wetland-nc",
        default=None,
        help="Path to peatland_halfdeg.nc (default: data/geodata_py/glwd3/)",
    )
    parser.add_argument(
        "--wetland-product",
        choices=["wforests", "noforests"],
        default="wforests",
        help=(
            "Which wetland product to use: 'wforests' includes forest wetland "
            "types (classes 5,6); 'noforests' excludes them (default: wforests)"
        ),
    )
    parser.add_argument(
        "--output-path",
        default=None,
        help="Output path (default: alongside LU file with _peatland tag)",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    from landsymm.config import get_remap_output_dir, get_geodata_dir

    lu_path = args.lu_path or str(
        get_remap_output_dir()
        / "remaps_v10_old_62892_gL" / "LU.remapv10_old_62892_gL.txt"
    )
    wetland_nc = args.wetland_nc or str(
        get_geodata_dir() / "glwd3" / "peatland_halfdeg.nc"
    )

    if args.output_path:
        output_path = args.output_path
    else:
        lu_dir = os.path.dirname(lu_path)
        lu_base = os.path.basename(lu_path)
        name, ext = os.path.splitext(lu_base)
        output_path = os.path.join(lu_dir, f"{name}_peatland{ext}")

    insert_wetland_approach_h(
        lu_path, wetland_nc, output_path,
        wetland_product=args.wetland_product,
        verbose=not args.quiet,
    )


if __name__ == "__main__":
    main()
