"""Insert GLWD3 wetland/peatland into scenario forLPJG landcover.txt files.

Applies the same Approach H logic as wetland_into_hilda.py, but operates on
the harmonized PLUMharm2LPJG output landcover.txt files found in each
scenario's forLPJG directory.

For each scenario, produces a new file 'landcover_peatland.txt' alongside the
original 'landcover.txt', with a PEATLAND column carved from NATURAL.

Approach H (recap):
  1. For each gridcell (Lon, Lat), find the minimum NATURAL fraction across
     ALL years in the dataset.
  2. PEATLAND = min(min_NATURAL, GLWD3_wetland_frac).
  3. NATURAL_new = NATURAL_old - PEATLAND at every year.
  4. CROPLAND, PASTURE, BARREN are unchanged.
  5. Fractions sum to 1.0 at every gridcell and year.
"""
from __future__ import annotations

import argparse
import glob
import os
import time

import numpy as np
import pandas as pd

from .wetland_into_hilda import (
    _load_wetland_frac_from_nc,
    _extract_wetland_frac_for_coords,
)


def _find_scenario_dirs(parent_dir: str) -> list[str]:
    """Discover scenario forLPJG directories containing landcover.txt."""
    pattern = os.path.join(parent_dir, "*", "s1.*.forLPJG")
    dirs = sorted(glob.glob(pattern))
    return [d for d in dirs if os.path.isfile(os.path.join(d, "landcover.txt"))]


def insert_peatland_into_landcover(
    landcover_path: str,
    wetland_nc_path: str,
    output_path: str,
    wetland_product: str = "wforests",
    verbose: bool = True,
) -> str:
    """Insert GLWD3 peatland into a forLPJG landcover.txt file.

    Parameters
    ----------
    landcover_path : str
        Path to the scenario's landcover.txt file.
    wetland_nc_path : str
        Path to the peatland_halfdeg.nc file.
    output_path : str
        Path for the output landcover_peatland.txt file.
    wetland_product : str
        "wforests" or "noforests".
    verbose : bool
        Print progress messages.

    Returns
    -------
    str
        Path to the output file.
    """
    t0 = time.perf_counter()

    if verbose:
        print(f"  Loading wetland fraction grid: {wetland_nc_path}")
        print(f"    Wetland product: {wetland_product}")
    wf_2d, lats_grid, lons_grid = _load_wetland_frac_from_nc(
        wetland_nc_path, wetland_product=wetland_product
    )

    if verbose:
        print(f"  Reading landcover: {landcover_path}")
    lc = pd.read_csv(landcover_path, sep=r"\s+")

    expected_cols = {"Lon", "Lat", "Year", "NATURAL", "CROPLAND", "PASTURE", "BARREN"}
    if not expected_cols.issubset(set(lc.columns)):
        missing = expected_cols - set(lc.columns)
        raise RuntimeError(f"Missing columns in landcover file: {missing}")

    n_rows = len(lc)
    years = sorted(lc["Year"].unique())
    n_years = len(years)
    if verbose:
        print(f"    {n_rows:,} rows, {n_years} years ({years[0]}-{years[-1]})")

    # Extract GLWD3 wetland fraction for each coordinate
    peatland_frac_raw = _extract_wetland_frac_for_coords(
        wf_2d, lats_grid, lons_grid,
        lc["Lon"].values, lc["Lat"].values,
    )

    # Approach H: minimum NATURAL per gridcell across all years
    if verbose:
        print("  Computing minimum NATURAL per gridcell...")
    min_nat = lc.groupby(["Lon", "Lat"])["NATURAL"].transform("min")

    peatland = np.minimum(min_nat.values, peatland_frac_raw)
    peatland = np.maximum(peatland, 0.0)

    new_natural = lc["NATURAL"].values - peatland

    n_neg = np.sum(new_natural < 0)
    if n_neg > 0:
        if verbose:
            print(f"  WARNING: {n_neg} negative NATURAL values, clamping to 0")
        neg_mask = new_natural < 0
        peatland[neg_mask] = lc["NATURAL"].values[neg_mask]
        new_natural[neg_mask] = 0.0

    # Build output — preserve original column order, insert PEATLAND before BARREN
    lc_out = lc.copy()
    lc_out["NATURAL"] = new_natural
    lc_out["PEATLAND"] = peatland

    orig_cols = list(lc.columns)
    barren_idx = orig_cols.index("BARREN")
    col_order = orig_cols[:barren_idx] + ["PEATLAND"] + orig_cols[barren_idx:]
    lc_out = lc_out[col_order]

    # Verify fractions sum to 1.0
    frac_cols = ["NATURAL", "CROPLAND", "PASTURE", "PEATLAND", "BARREN"]
    row_sums = lc_out[frac_cols].sum(axis=1)
    max_dev = np.max(np.abs(row_sums - 1.0))
    if verbose:
        print(f"  Max deviation from sum=1.0: {max_dev:.2e}")
    if max_dev > 1e-4:
        print(f"  WARNING: Row sums deviate from 1.0 by up to {max_dev:.6e}")

    # Summary
    if verbose:
        n_with = np.sum(peatland > 0)
        n_cells = len(lc_out.loc[lc_out["PEATLAND"] > 0, ["Lon", "Lat"]].drop_duplicates())
        print(f"  Rows with PEATLAND > 0: {n_with:,} / {n_rows:,}")
        print(f"  Unique gridcells with PEATLAND > 0: {n_cells:,}")
        if n_with > 0:
            print(f"  PEATLAND range: [{peatland[peatland > 0].min():.6f}, {peatland.max():.6f}]")

    # Write output — vectorized for speed on large files (~5M rows)
    if verbose:
        print(f"  Writing: {output_path}")

    header_line = " ".join(col_order)
    fmt_parts = []
    for c in col_order:
        if c in ("Lon", "Lat"):
            fmt_parts.append("%.2f")
        elif c == "Year":
            fmt_parts.append("%d")
        else:
            fmt_parts.append("%.6f")

    np.savetxt(
        output_path,
        lc_out.values,
        fmt=fmt_parts,
        delimiter=" ",
        header=header_line,
        comments="",
    )

    elapsed = time.perf_counter() - t0
    if verbose:
        print(f"  Done ({elapsed:.1f}s)")

    return output_path


def main(
    parent_dir: str | None = None,
    wetland_nc_path: str | None = None,
    wetland_product: str = "wforests",
    scenarios: list[str] | None = None,
    verbose: bool = True,
) -> None:
    """Process all scenario forLPJG landcover.txt files."""
    from landsymm.config import get_plum_output_dir, get_geodata_dir

    t0_total = time.perf_counter()

    if parent_dir is None:
        parent_dir = str(get_plum_output_dir())
    if wetland_nc_path is None:
        wetland_nc_path = str(
            get_geodata_dir() / "glwd3" / "peatland_halfdeg.nc"
        )

    if not os.path.isfile(wetland_nc_path):
        raise FileNotFoundError(
            f"Wetland NetCDF not found: {wetland_nc_path}\n"
            "Run glwd3_to_halfdeg.py first to generate peatland_halfdeg.nc"
        )

    all_dirs = _find_scenario_dirs(parent_dir)
    if scenarios:
        all_dirs = [d for d in all_dirs if any(s in d for s in scenarios)]

    if not all_dirs:
        print(f"No forLPJG directories found under {parent_dir}")
        return

    print("=" * 60)
    print("  Insert PEATLAND into scenario forLPJG landcover.txt files")
    print(f"  Parent directory: {parent_dir}")
    print(f"  Wetland product: {wetland_product}")
    print(f"  Found {len(all_dirs)} scenario(s)")
    print("=" * 60)

    for i, fdir in enumerate(all_dirs, 1):
        scenario = os.path.basename(os.path.dirname(fdir))
        lc_path = os.path.join(fdir, "landcover.txt")
        out_path = os.path.join(fdir, "landcover_peatland.txt")

        print(f"\n[{i}/{len(all_dirs)}] {scenario}")
        insert_peatland_into_landcover(
            lc_path, wetland_nc_path, out_path,
            wetland_product=wetland_product,
            verbose=verbose,
        )

    elapsed = time.perf_counter() - t0_total
    minutes = int(elapsed // 60)
    seconds = int(elapsed % 60)
    print(f"\n{'=' * 60}")
    print(f"  All scenarios processed ({minutes}m {seconds}s)")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Insert GLWD3 peatland into scenario forLPJG landcover.txt files. "
            "Produces landcover_peatland.txt alongside each original."
        )
    )
    parser.add_argument(
        "--parent-dir", default=None,
        help="Parent directory containing scenario folders (default: data/PLUMv2_LU_default_output)",
    )
    parser.add_argument(
        "--wetland-nc", default=None,
        help="Path to peatland_halfdeg.nc (default: data/geodata_py/glwd3/peatland_halfdeg.nc)",
    )
    parser.add_argument(
        "--wetland-product",
        choices=["wforests", "noforests"],
        default="wforests",
        help="Which wetland product to use (default: wforests)",
    )
    parser.add_argument(
        "--scenarios", nargs="*", default=None,
        help="Specific scenarios to process (default: all found)",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    main(
        parent_dir=args.parent_dir,
        wetland_nc_path=args.wetland_nc,
        wetland_product=args.wetland_product,
        scenarios=args.scenarios,
        verbose=not args.quiet,
    )
