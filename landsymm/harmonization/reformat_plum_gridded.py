"""Reformat raw PLUM gridded outputs into harmonization-ready inputs.

Python port of reformat_gridded_updated.R.

Converts raw PLUM LandCover.txt and LandUse.txt files into the format
expected by PLUMharm:
  - LandCoverFract.txt  (NATURAL, CROPLAND, PASTURE, BARREN, URBAN fractions)
  - CropFract.txt       (per-crop fractions of cropland)
  - Fert.txt            (per-crop fertilization rates)
  - Irrig.txt           (per-crop irrigation intensity)
  - LandUse.txt         (reformatted combined table)

The original LandUse.txt is renamed to LandUse_orig.txt(.gz) before
the reformatted version is written.
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

CROP_ORDER = [
    "pasture", "energycrops", "wheat", "oilcropsNFix", "oilcropsOther",
    "maize", "rice", "pulses", "setaside", "starchyRoots", "fruitveg", "sugar",
]

VAR_ORDER = ["A", "FI", "FQ", "II", "IQ", "OI", "Y"]

CROP_FRACT_RENAME = {
    "wheat_A": "CerealsC3",
    "oilcropsNFix_A": "OilcropsNFix",
    "oilcropsOther_A": "OilcropsOther",
    "starchyRoots_A": "StarchyRoots",
    "pulses_A": "Pulses",
    "maize_A": "Maize",
    "energycrops_A": "Miscanthus",
    "rice_A": "Rice",
    "fruitveg_A": "FruitVeg",
    "sugar_A": "Sugar",
    "setaside_A": "Setaside",
}

FERT_RENAME = {k.replace("_A", "_FQ"): v for k, v in CROP_FRACT_RENAME.items()}
IRRIG_RENAME = {k.replace("_A", "_II"): v for k, v in CROP_FRACT_RENAME.items()}


def _find_file(data_dir: str, basename: str) -> str:
    gz_path = os.path.join(data_dir, f"{basename}.gz")
    txt_path = os.path.join(data_dir, basename)
    if os.path.exists(gz_path):
        return gz_path
    if os.path.exists(txt_path):
        return txt_path
    raise FileNotFoundError(f"Neither {basename} nor {basename}.gz found in {data_dir}")


def _read_table(path: str) -> pd.DataFrame:
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return pd.read_csv(f, sep=r"\s+")
    return pd.read_csv(path, sep=r"\s+")


def _format_df(df: pd.DataFrame, lon_lat_prec: int = 2, data_prec: int = 10) -> pd.DataFrame:
    out = df.copy()
    for col in out.columns:
        if col in ("Lon", "Lat"):
            out[col] = out[col].map(lambda x: f"{x:.{lon_lat_prec}f}")
        else:
            out[col] = out[col].map(lambda x: f"{x:.{data_prec}f}")
    return out


def _write_table(df: pd.DataFrame, path: str, lon_lat_prec: int = 2, data_prec: int = 10) -> None:
    formatted = _format_df(df, lon_lat_prec, data_prec)
    formatted.to_csv(path, sep=" ", index=False, quoting=3)


def reformat_gridded(data_dir: str, verbose: bool = True) -> None:
    """Reformat a single year's PLUM outputs."""
    if verbose:
        print(f"  Reformatting {data_dir}")

    year = int(os.path.basename(data_dir))

    lu_fpath = _find_file(data_dir, "LandUse.txt")
    lc_fpath = _find_file(data_dir, "LandCover.txt")

    # --- Process LandCover ---
    lc_raw = _read_table(lc_fpath)

    prot_df = lc_raw[["Lon", "Lat", "Protection", "TotalArea"]].pivot_table(
        index=["Lon", "Lat"], columns="Protection", values="TotalArea", aggfunc="sum"
    ).reset_index()
    for col in ["protected", "unprotected"]:
        if col not in prot_df.columns:
            prot_df[col] = 0.0
    prot_df["total"] = prot_df["protected"] + prot_df["unprotected"]
    prot_df["pa_fraction"] = np.where(prot_df["total"] > 0, prot_df["protected"] / prot_df["total"], 0.0)
    prot_df = prot_df[["Lon", "Lat", "protected", "pa_fraction"]]

    lc_agg = lc_raw.groupby(["Lon", "Lat"], as_index=False).agg(
        {c: "sum" for c in lc_raw.columns if c not in ("Lon", "Lat", "Protection")}
    )
    lc_agg["suitable"] = (
        lc_agg["Cropland"] + lc_agg["Pasture"] + lc_agg["TimberForest"]
        + lc_agg["UnmanagedForest"] + lc_agg["CarbonForest"]
        + lc_agg["OtherNatural"] + lc_agg.get("Photovoltaics", 0) + lc_agg.get("Agrivoltaics", 0)
    )
    lc_agg["otherNatural"] = lc_agg["OtherNatural"] + lc_agg.get("Photovoltaics", 0)
    lc_agg["cropland"] = lc_agg["Cropland"] + lc_agg.get("Agrivoltaics", 0)

    lc_df = lc_agg.merge(prot_df, on=["Lon", "Lat"], how="left")
    lc_df = lc_df.rename(columns={
        "TotalArea": "area",
        "TimberForest": "timberForest",
        "CarbonForest": "carbonForest",
        "UnmanagedForest": "unmanagedForest",
        "Pasture": "pasture",
        "Barren": "barren",
        "Urban": "urban",
    })
    keep_cols = [
        "Lon", "Lat", "area", "suitable", "protected", "pa_fraction",
        "timberForest", "carbonForest", "unmanagedForest", "otherNatural",
        "cropland", "pasture", "barren", "urban",
    ]
    lc_df = lc_df[[c for c in keep_cols if c in lc_df.columns]]

    # --- Process LandUse ---
    lu_raw = _read_table(lu_fpath)

    mgmt_cols = [c for c in ["FI", "FQ", "II", "IQ", "OI", "Y"] if c in lu_raw.columns]
    for col in mgmt_cols:
        lu_raw[col] = lu_raw[col] * lu_raw["A"]

    crop_df = lu_raw.groupby(["Lon", "Lat", "Crop"], as_index=False).agg(
        {c: "sum" for c in ["A"] + mgmt_cols}
    )
    for col in mgmt_cols:
        crop_df[col] = np.where(crop_df["A"] > 0, crop_df[col] / crop_df["A"], 0.0)

    crop_df = crop_df.merge(lc_df[["Lon", "Lat", "cropland"]], on=["Lon", "Lat"], how="left")
    crop_df["A"] = np.where(
        crop_df["Crop"] == "pasture", 1.0,
        np.where(crop_df["cropland"] > 0, crop_df["A"] / crop_df["cropland"], 0.0)
    )
    crop_df = crop_df.drop(columns=["cropland"])

    melted = crop_df.melt(id_vars=["Lon", "Lat", "Crop"], value_vars=["A"] + mgmt_cols)
    melted["Var"] = melted["Crop"] + "_" + melted["variable"]
    crop_order_map = {c: i for i, c in enumerate(CROP_ORDER)}
    var_order_map = {v: i for i, v in enumerate(VAR_ORDER)}
    melted["crop_rank"] = melted["Crop"].map(crop_order_map).fillna(999)
    melted["var_rank"] = melted["variable"].map(var_order_map).fillna(999)
    melted = melted.sort_values(["crop_rank", "var_rank"])
    melted = melted.drop(columns=["variable", "Crop", "crop_rank", "var_rank"])
    crop_wide = melted.pivot_table(index=["Lon", "Lat"], columns="Var", values="value", aggfunc="first").reset_index()
    crop_wide.columns.name = None

    # Rename original LandUse file
    if lu_fpath.endswith(".gz"):
        new_name = os.path.join(data_dir, "LandUse_orig.txt.gz")
    else:
        new_name = os.path.join(data_dir, "LandUse_orig.txt")
    if not os.path.exists(new_name):
        os.rename(lu_fpath, new_name)

    # Join land cover and crop data
    lu_df = lc_df.merge(crop_wide, on=["Lon", "Lat"], how="left")
    lu_df = lu_df[lu_df["area"] > 0].copy()
    lu_df = lu_df.fillna(0.0)
    lu_df["managedForest"] = lu_df["timberForest"] + lu_df["carbonForest"]

    # --- Write CropFract.txt ---
    a_cols = {k: v for k, v in CROP_FRACT_RENAME.items() if k in lu_df.columns}
    cf_df = lu_df[["Lon", "Lat"]].copy()
    for src, dst in a_cols.items():
        cf_df[dst] = lu_df[src]
    _write_table(cf_df, os.path.join(data_dir, "CropFract.txt"))

    # --- Write Fert.txt ---
    fq_cols = {k: v for k, v in FERT_RENAME.items() if k in lu_df.columns}
    fert_df = lu_df[["Lon", "Lat"]].copy()
    for src, dst in fq_cols.items():
        fert_df[dst] = lu_df[src]
    _write_table(fert_df, os.path.join(data_dir, "Fert.txt"))

    # --- Write Irrig.txt ---
    ii_cols = {k: v for k, v in IRRIG_RENAME.items() if k in lu_df.columns}
    irrig_df = lu_df[["Lon", "Lat"]].copy()
    for src, dst in ii_cols.items():
        irrig_df[dst] = lu_df[src]
    _write_table(irrig_df, os.path.join(data_dir, "Irrig.txt"))

    # --- Write LandCoverFract.txt ---
    lcf_df = lu_df[["Lon", "Lat"]].copy()
    lcf_df["NATURAL"] = (lu_df["managedForest"] + lu_df["unmanagedForest"] + lu_df["otherNatural"]) / lu_df["area"]
    lcf_df["CROPLAND"] = lu_df["cropland"] / lu_df["area"]
    lcf_df["PASTURE"] = lu_df["pasture"] / lu_df["area"]
    lcf_df["BARREN"] = lu_df["barren"] / lu_df["area"]
    lcf_df["URBAN"] = lu_df["urban"] / lu_df["area"]
    _write_table(lcf_df, os.path.join(data_dir, "LandCoverFract.txt"))

    # --- Write reformatted LandUse.txt ---
    lu_out_cols = ["Lon", "Lat", "area", "suitable", "protected", "pa_fraction",
                   "managedForest", "unmanagedForest", "otherNatural", "cropland",
                   "pasture", "barren", "urban"]
    for crop in CROP_ORDER:
        for var in VAR_ORDER:
            col = f"{crop}_{var}"
            if col in lu_df.columns:
                lu_out_cols.append(col)
    lu_out = lu_df[[c for c in lu_out_cols if c in lu_df.columns]]
    _write_table(lu_out, os.path.join(data_dir, "LandUse.txt"))


def reformat_scenario(scenario_dir: str, verbose: bool = True) -> None:
    """Reformat all years for a single scenario (e.g., SSP1_RCP26)."""
    s1_dir = os.path.join(scenario_dir, "s1")
    if not os.path.isdir(s1_dir):
        raise RuntimeError(f"Directory {s1_dir} does not exist")

    year_dirs = sorted([
        d for d in os.listdir(s1_dir)
        if os.path.isdir(os.path.join(s1_dir, d)) and re.match(r"^\d{4}$", d)
    ])

    if verbose:
        print(f"Processing {scenario_dir}: {len(year_dirs)} years")

    for yr_dir_name in year_dirs:
        reformat_gridded(os.path.join(s1_dir, yr_dir_name), verbose=verbose)


def reformat_all(
    parent_dir: str,
    scenarios: Sequence[str] | None = None,
    verbose: bool = True,
) -> None:
    """Reformat all scenarios under the parent directory."""
    if scenarios is None:
        all_dirs = sorted([
            d for d in os.listdir(parent_dir)
            if os.path.isdir(os.path.join(parent_dir, d))
            and re.match(r"^SSP[1-5]_RCP\d{2}$", d)
        ])
    else:
        all_dirs = list(scenarios)

    for ssp in all_dirs:
        scenario_dir = os.path.join(parent_dir, ssp)
        if not os.path.isdir(scenario_dir):
            print(f"Skipping {ssp}: directory not found")
            continue
        reformat_scenario(scenario_dir, verbose=verbose)

    if verbose:
        print("All scenarios reformatted.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reformat raw PLUM gridded outputs for PLUMharm (Python port of reformat_gridded_updated.R)."
    )
    parser.add_argument(
        "--parent-dir",
        default=None,
        help="Parent directory containing SSP scenario dirs (default: data/PLUMv2_LU_default_output)",
    )
    parser.add_argument(
        "--scenarios", nargs="+", default=None,
        help="Scenarios to process (default: all matching SSP*_RCP*)",
    )
    parser.add_argument(
        "--single-year-dir",
        default=None,
        help="Process a single year directory (e.g., .../SSP1_RCP26/SSP1_RCP26/s1/2021)",
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress progress output")
    args = parser.parse_args()

    if args.single_year_dir:
        reformat_gridded(args.single_year_dir, verbose=not args.quiet)
    else:
        if args.parent_dir is None:
            from landsymm.config import get_plum_output_dir
            parent_dir = str(get_plum_output_dir())
        else:
            parent_dir = args.parent_dir
        reformat_all(parent_dir, scenarios=args.scenarios, verbose=not args.quiet)


if __name__ == "__main__":
    main()
