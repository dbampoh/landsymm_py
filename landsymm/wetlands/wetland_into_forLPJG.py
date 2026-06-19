"""Insert GLWD3 wetland/peatland into scenario forLPJG landcover.txt files.

Why this module exists (OPTIONAL stage)
=======================================
This is the scenario-side companion of ``wetland_into_hilda.py`` and
belongs to **Stage 4 of the landsymm_py pipeline, which is optional**.
Run it only if your downstream LPJ-GUESS runs need an explicit
PEATLAND land-cover class — typically because you are running coupled
LPJ-GUESS ↔ IMOGEN climate simulations where peatland CH₄ emissions
feed back into the IMOGEN climate model alongside CO₂ and N₂O.

If your downstream LPJ-GUESS runs use prescribed (offline) climate
forcing and do not need explicit peatland CH₄ accounting, you can
skip this stage — the Stage-3 PLUMharm2LPJG output landcover.txt
files are already complete LPJ-GUESS-ready scenario inputs.

What this module does
=====================
Applies the same Approach H logic as ``wetland_into_hilda.py`` but
operates on the harmonized PLUMharm2LPJG output ``landcover.txt`` files
found in each scenario's forLPJG directory. For each scenario, produces
a new file ``landcover_peatland.txt`` alongside the original
``landcover.txt``, with a PEATLAND column carved from NATURAL. The
``_peatland`` suffix keeps peatland and non-peatland scenario LU files
separate so a single project can support both kinds of LPJ-GUESS runs
side by side.

The min-over-years approach guarantees that PEATLAND is bounded by the
minimum natural-vegetation availability across all years in the
scenario, so it never exceeds available NATURAL fraction at any year.

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


def _find_scenario_dirs(parent_dir: str, member: str | None = None) -> list[str]:
    """Discover scenario forLPJG directories containing landcover.txt."""
    from landsymm.config import get_member

    if member is None:
        member = get_member()
    pattern = os.path.join(parent_dir, "*", f"{member}.*.forLPJG")
    dirs = sorted(glob.glob(pattern))
    return [d for d in dirs if os.path.isfile(os.path.join(d, "landcover.txt"))]


def _resolve_landcover_files(
    parent_dir: str | None = None,
    flat_parent: str | None = None,
    landcover: list[str] | None = None,
    scenarios: list[str] | None = None,
) -> list[str]:
    """Resolve the list of ``landcover.txt`` files to process (first match wins):

    - ``landcover``: explicit ``landcover.txt`` file paths, used as-is.
    - ``flat_parent``: scan ``<flat_parent>/*/landcover.txt`` (flat layout, e.g.
      ``Data/lu/<scenario>/landcover.txt``).
    - ``parent_dir``: default forLPJG scan,
      ``<parent_dir>/*/<member>.*.forLPJG/landcover.txt``.

    The ``scenarios`` substring filter applies to the scan modes only (not to an
    explicit ``landcover`` list).
    """
    if landcover:
        files = [os.path.abspath(p) for p in landcover]
        missing = [p for p in files if not os.path.isfile(p)]
        if missing:
            raise FileNotFoundError(
                "landcover file(s) not found: " + ", ".join(missing)
            )
        return files
    if flat_parent:
        files = sorted(glob.glob(os.path.join(flat_parent, "*", "landcover.txt")))
    else:
        files = [
            os.path.join(d, "landcover.txt") for d in _find_scenario_dirs(parent_dir)
        ]
    if scenarios:
        files = [f for f in files if any(s in f for s in scenarios)]
    return files


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
    landcover: list[str] | None = None,
    flat_parent: str | None = None,
) -> None:
    """Insert PEATLAND into scenario landcover.txt file(s).

    By default scans the PLUM parent dir for forLPJG scenario directories. Use
    ``landcover`` to pass explicit landcover.txt paths, or ``flat_parent`` to scan
    a flat ``<dir>/<scenario>/landcover.txt`` layout (e.g. a paper ``Data/lu/``).
    A ``landcover_peatland.txt`` is written alongside each input.
    """
    from landsymm.config import get_geodata_dir, get_plum_output_dir

    t0_total = time.perf_counter()

    if wetland_nc_path is None:
        wetland_nc_path = str(get_geodata_dir() / "glwd3" / "peatland_halfdeg.nc")
    if not os.path.isfile(wetland_nc_path):
        raise FileNotFoundError(
            f"Wetland NetCDF not found: {wetland_nc_path}\n"
            "Run glwd3_to_halfdeg.py first to generate peatland_halfdeg.nc"
        )

    # Default to the PLUM parent dir only when no explicit files / flat parent given.
    if landcover is None and flat_parent is None and parent_dir is None:
        parent_dir = str(get_plum_output_dir())

    if landcover:
        mode = f"explicit ({len(landcover)} file(s))"
    elif flat_parent:
        mode = f"flat scan: {flat_parent}/*/landcover.txt"
    else:
        mode = f"forLPJG scan: {parent_dir}"

    files = _resolve_landcover_files(
        parent_dir=parent_dir, flat_parent=flat_parent,
        landcover=landcover, scenarios=scenarios,
    )
    if not files:
        print(f"No landcover.txt files found ({mode})")
        return

    print("=" * 60)
    print("  Insert PEATLAND into scenario landcover.txt files")
    print(f"  Mode: {mode}")
    print(f"  Wetland product: {wetland_product}")
    print(f"  Found {len(files)} file(s)")
    print("=" * 60)

    for i, lc_path in enumerate(files, 1):
        out_path = os.path.join(os.path.dirname(lc_path), "landcover_peatland.txt")
        label = os.path.basename(os.path.dirname(lc_path)) or lc_path

        print(f"\n[{i}/{len(files)}] {label}")
        insert_peatland_into_landcover(
            lc_path, wetland_nc_path, out_path,
            wetland_product=wetland_product,
            verbose=verbose,
        )

    elapsed = time.perf_counter() - t0_total
    minutes = int(elapsed // 60)
    seconds = int(elapsed % 60)
    print(f"\n{'=' * 60}")
    print(f"  All files processed ({minutes}m {seconds}s)")
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
        help="Parent dir for the default forLPJG scan "
             "(<parent>/*/<member>.*.forLPJG/landcover.txt; default: config.get_plum_output_dir()).",
    )
    parser.add_argument(
        "--flat-parent", default=None,
        help="Scan a flat layout instead: <flat-parent>/*/landcover.txt "
             "(e.g. a paper Data/lu/ directory).",
    )
    parser.add_argument(
        "--landcover", nargs="+", default=None,
        help="Explicit landcover.txt file path(s) to process (overrides the scans).",
    )
    parser.add_argument(
        "--wetland-nc", default=None,
        help="Path to peatland_halfdeg.nc (default: <geodata>/glwd3/peatland_halfdeg.nc)",
    )
    parser.add_argument(
        "--wetland-product",
        choices=["wforests", "noforests"],
        default="wforests",
        help="Which wetland product to use (default: wforests)",
    )
    parser.add_argument(
        "--scenarios", nargs="*", default=None,
        help="Substring filter for the scan modes (ignored with --landcover).",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    main(
        parent_dir=args.parent_dir,
        wetland_nc_path=args.wetland_nc,
        wetland_product=args.wetland_product,
        scenarios=args.scenarios,
        verbose=not args.quiet,
        landcover=args.landcover,
        flat_parent=args.flat_parent,
    )
