"""Main driver for remapping pipeline (stub).

Pseudocode:
1) Load config (remap_options).
2) Setup output dirs + logging.
3) Read output gridlist.
4) Interpolate soil inputs.
5) Optional climate gridlist diagnostics.
6) Import LU (HILDA+ or LUH2), handle unveg, normalize to 1.
7) Import MIRCA crop fractions + mapping (get_remapv2_keys).
8) Import Nfert (AgGRID) and align with crops.
9) Interpolate cropfracs/nfert to gridlist; normalize.
10) Write LU/cropfracs/nfert/gridlist outputs.
"""

# TODO:
# - implement full remap pipeline with exact parity
# - wire logging/diary to file
# - ensure output formatting matches MATLAB
from __future__ import annotations

import atexit
import os
import sys
import warnings
from dataclasses import dataclass
from typing import TextIO

# Avoid Qt/GUI warnings during non-interactive runs
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
from netCDF4 import Dataset

from landsymm.common.interpolation import inpaint_nans
from landsymm.common.io_utils import delete_existing_outputs
from landsymm.common.lpjg_io import save_table
from landsymm.common.mapping_tools import xz_to_yxz

from .cropfrac import process_mirca
from .diagnostics import map_missing
from .grid import load_gridlist
from .lu_import import import_lu_hildaplus, import_lu_luh2, note_unveg
from .nfert import aggrid_type_to_landsymm, combine_nfert
from .remap_options import RemapConfig
from .soil import interpolate_soil


class RemapPipeline:
    """Class-based remapping pipeline (skeleton)."""

    def __init__(self, cfg: RemapConfig) -> None:
        self.cfg = cfg

    def run(self) -> None:
        """Run the remapping workflow (TODO)."""
        run_remap(self.cfg)


def run_remap(cfg: RemapConfig) -> None:
    """Run remapping pipeline using config (TODO)."""
    warnings.filterwarnings(
        "ignore",
        message="force_mat_save and force_mat_nosave can't both be true.*",
        category=UserWarning,
    )
    warnings.filterwarnings(
        "ignore",
        message="MAT table missing Lon/Lat headers; reading text instead.",
        category=UserWarning,
    )
    warnings.filterwarnings(
        "ignore",
        message="Assuming columns Lon, Lat.",
        category=UserWarning,
    )
    warnings.filterwarnings(
        "ignore",
        message="invalid value encountered in divide",
        category=RuntimeWarning,
    )

    out_dir = os.path.join(cfg.out_dir_top, f"remaps_v{cfg.remap_ver}")
    os.makedirs(out_dir, exist_ok=True)
    out_dir_figs = os.path.join(out_dir, "figs")
    os.makedirs(out_dir_figs, exist_ok=True)

    log_path = cfg.log_file or os.path.join(out_dir, "remap.log")
    log_file = open(log_path, "w")
    stdout_prev = sys.stdout
    stderr_prev = sys.stderr
    sys.stdout = _TeeWriter(stdout_prev, log_file)
    sys.stderr = _TeeWriter(stderr_prev, log_file)

    def _restore_streams() -> None:
        sys.stdout = stdout_prev
        sys.stderr = stderr_prev
        if not log_file.closed:
            log_file.close()

    atexit.register(_restore_streams)

    print(f"Log file: {log_path}")
    print(f"PLUMsetAside_frac: {cfg.plum_setaside_frac}")
    print(f"inpaint_method: {cfg.inpaint_method}")
    print(f"year_list_out: {min(cfg.year_list_out)}-{max(cfg.year_list_out)}")
    print(f"this_ver: {cfg.this_ver}")
    print(f"remap_ver: {cfg.remap_ver}")
    print(f"force_all_rainfed: {int(cfg.force_all_rainfed)}")

    n_years_out = len(cfg.year_list_out)
    interp_test_y1 = cfg.interp_test_y1 or cfg.year_list_out[0]
    interp_test_yN = cfg.interp_test_yN or cfg.year_list_out[-1]
    year_list_interp_test = list(range(interp_test_y1, interp_test_yN + 1))
    year_list_out = list(cfg.year_list_out)
    if not set(year_list_interp_test).intersection(year_list_out):
        print(
            f"Warning: interpolation test years ({interp_test_y1}-{interp_test_yN}) outside year_list_out "
            f"({year_list_out[0]}-{year_list_out[-1]}). Using entire year_list_out."
        )
        interp_test_y1 = year_list_out[0]
        interp_test_yN = year_list_out[-1]
        year_list_interp_test = year_list_out
    I_yrs_in_test_period = [year_list_out.index(y) for y in year_list_interp_test if y in year_list_out]
    interp_test_period_str = f"{interp_test_y1}-{interp_test_yN}"

    # Output gridlist
    gridlist = load_gridlist(cfg.file_gridlist_out)
    n_cells = len(gridlist["list2map"])

    # Soil interpolation
    interpolate_soil(
        cfg.files_soil,
        gridlist,
        out_dir,
        out_dir_figs,
        cfg.inpaint_method,
        cfg.remap_ver,
        cfg.out_width,
        cfg.delimiter,
        cfg.overwrite,
        cfg.fancy,
    )

    # Climate gridlist diagnostics
    if cfg.file_gridlist_climate and os.path.exists(cfg.file_gridlist_climate):
        climate_gridlist = load_gridlist(cfg.file_gridlist_climate)
        missing = np.setdiff1d(gridlist["list2map"], climate_gridlist["list2map"])
        if missing.size > 0:
            missing_lonlats = gridlist["lonlats"][np.isin(gridlist["list2map"], missing)]
            np.savetxt(os.path.join(out_dir, "missing_climate.csv"), missing_lonlats, delimiter=",")

    # Import land uses
    if cfg.lu_source == "LUH2":
        out_lu, carea_yx = import_lu_luh2(cfg.geodata_dir, cfg.year_list_out, gridlist)
    elif cfg.lu_source == "HILDA+":
        out_lu, carea_yx = import_lu_hildaplus(
            cfg.geodata_dir, cfg.year_list_out, gridlist,
            max_nan_frac=cfg.max_nan_frac,
        )
    else:
        raise RuntimeError(f"lu_source {cfg.lu_source} not recognized")

    # Note unvegetated cells
    out_lu = note_unveg(out_lu, cfg.fill_unveg)

    # Force LU sum to 1
    lu_out_x1y = np.sum(out_lu["garr_xvy"], axis=1)
    j = 0
    n_lu_out = len(out_lu["varNames"])
    while _max_sum_diff(out_lu["garr_xvy"], 1.0) > 1e-6:
        if j == 0:
            print("Forcing all land cells' LU fractions to sum to 1:")
            print(f"    Before iteration 1, max err = {np.max(np.abs(lu_out_x1y - 1))}")
        j += 1
        if j > 50:
            raise RuntimeError('Possible infinite loop in "remap_import_lu_force_sum1()".')
        out_lu["garr_xvy"] = out_lu["garr_xvy"] / np.tile(lu_out_x1y.reshape(n_cells, 1, n_years_out), (1, n_lu_out, 1))
        lu_out_x1y = np.sum(out_lu["garr_xvy"], axis=1)
        print(f"    After  iteration {j}, max err = {np.max(np.abs(lu_out_x1y - 1))}")
    if cfg.fill_unveg > 0:
        note_unveg(out_lu, -1)

    # Crop fractions
    file_cropmirca = os.path.join(
        cfg.geodata_dir, "MIRCA", "harvested_area_grids_26crops_30mn", "MIRCA.txt"
    )
    (
        mid_cropfrac,
        croparea_in,
        list_crops_out,
        list_crops_combined_out,
        list_ignore_frac,
        in2out_key_frac,
        all_ver_names,
        all_ver_ignore_types,
    ) = process_mirca(
        file_cropmirca,
        gridlist,
        cfg.this_ver,
        cfg.force_all_rainfed,
        cfg.plum_setaside_frac,
    )

    ncrops_out = len(list_crops_out)
    list_crops_out_nfert = [c for c in list_crops_out if c != "ExtraCrop"]

    # Nfert import
    print("Importing Nfert...")
    nfert_dir = os.path.join(cfg.geodata_dir, "AgGRID_nutrient_input_v1.1")
    print(f"nfert_dir: {nfert_dir}")

    list_crops_out_as_nfert = []
    for crop in list_crops_out_nfert:
        base = crop[:-1] if crop.endswith("i") else crop
        list_crops_out_as_nfert.append(aggrid_type_to_landsymm(base))

    mid_nfert = {"varNames": list_crops_out_nfert, "garr_xv": np.full((n_cells, len(list_crops_out_nfert)), np.nan)}
    list_crops_combined_nfert_in = sorted(set(list_crops_out_as_nfert))
    for crop in list_crops_combined_nfert_in:
        if "combined" in crop:
            in_x = combine_nfert(crop, croparea_in, nfert_dir, gridlist)
        else:
            this_file = os.path.join(nfert_dir, f"agmip_{crop}_apprate_fill_NPK_0.5.nc4")
            with Dataset(this_file, "r") as ds:
                var = ds.variables["Napprate"][:]
                if np.ma.isMaskedArray(var):
                    var = var.filled(np.nan)
                in_yx = np.flipud(np.asarray(var))
            in_x = in_yx.ravel(order="F")[np.asarray(gridlist["list2map"]) - 1]
        inds = [i for i, v in enumerate(list_crops_out_as_nfert) if v == crop]
        mid_nfert["garr_xv"][:, inds] = np.tile(in_x.reshape(-1, 1), (1, len(inds)))

    mid_nfert["list2map"] = gridlist["list2map"]
    mid_nfert["lonlats"] = gridlist["lonlats"]

    # Check RF/IR equality
    ir_inds = []
    for i, crop in enumerate(mid_nfert["varNames"]):
        if crop.endswith("i"):
            ir_inds.append(i)
    is_ir = np.zeros(len(mid_nfert["varNames"]), dtype=bool)
    is_ir[ir_inds] = True
    rf_inds = np.where(~is_ir)[0]
    for i in rf_inds:
        crop = mid_nfert["varNames"][i]
        crop_i = f"{crop}i"
        if crop_i in mid_nfert["varNames"]:
            j = mid_nfert["varNames"].index(crop_i)
            tmp_rf = mid_nfert["garr_xv"][:, i].copy()
            tmp_ir = mid_nfert["garr_xv"][:, j].copy()
            tmp_rf[np.isnan(tmp_rf)] = -1
            tmp_ir[np.isnan(tmp_ir)] = -1
            nbad = np.sum(tmp_rf != tmp_ir)
            if nbad > 0:
                raise RuntimeError(f"{crop} ({i},{j}): {nbad}")
    print("Rainfed and irrigated fertilization is equal.")

    # Convert kg/ha to kg/m2
    mid_nfert["garr_xv"] = mid_nfert["garr_xv"] * 1e-4
    print("Done.")

    # Interpolate cropfrac/nfert to gridlist
    print("Interpolating to match gridlist...")
    isbad_yx = map_missing(out_lu["garr_xvy"], gridlist, "LU", out_dir_figs)
    if np.any(isbad_yx[gridlist["mask_YX"]]):
        raise RuntimeError("You also need to interpolate land use")
    map_missing(mid_cropfrac["garr_xv"], gridlist, "cropfrac", out_dir_figs)
    map_missing(mid_nfert["garr_xv"], gridlist, "nfert", out_dir_figs)

    out_cropfrac = {**mid_cropfrac, "garr_xv": np.full((n_cells, ncrops_out), np.nan)}
    out_nfert = {**mid_cropfrac, "garr_xv": np.full((n_cells, ncrops_out), np.nan)}

    test_croparea_x = np.mean(
        out_lu["garr_xvy"][:, out_lu["varNames"].index("CROPLAND"), :][:, I_yrs_in_test_period],
        axis=1,
    ) * carea_yx.ravel(order="F")[np.asarray(gridlist["list2map"]) - 1]

    for c, crop in enumerate(list_crops_out):
        print(f"Interpolating {crop} ({c+1} of {ncrops_out})...")
        # cropfrac
        tmp0_yx = xz_to_yxz(mid_cropfrac["garr_xv"][:, c:c+1], gridlist["mask_YX"].shape, gridlist["list2map"])[:, :, 0]
        tmp1_yx = inpaint_nans(tmp0_yx, cfg.inpaint_method)
        tmp0_yx[np.isnan(tmp0_yx)] = 0
        area0 = np.sum(test_croparea_x * tmp0_yx.ravel(order="F")[np.asarray(gridlist["list2map"]) - 1])
        area1 = np.sum(test_croparea_x * tmp1_yx.ravel(order="F")[np.asarray(gridlist["list2map"]) - 1])
        area_diff = area1 - area0
        area_diff_pct = area_diff / area0 * 100 if area0 != 0 else 0
        print(f"    Area increased {area_diff_pct:0.2f}% ({area_diff*100:0.1g} ha/yr) in {interp_test_period_str}")
        out_cropfrac["garr_xv"][:, c] = tmp1_yx.ravel(order="F")[np.asarray(gridlist["list2map"]) - 1]
        if np.isnan(out_cropfrac["garr_xv"][:, c]).any():
            raise RuntimeError("NaN remaining in out_cropfrac.garr_xv(:,c)!")

        # nfert
        c2 = c
        c3 = c if crop in mid_nfert["varNames"] else None
        if c3 is not None:
            tmp0_yx = xz_to_yxz(mid_nfert["garr_xv"][:, c3:c3+1], gridlist["mask_YX"].shape, gridlist["list2map"])[:, :, 0]
        else:
            tmp0_yx = np.zeros(gridlist["mask_YX"].shape)
        tmp1_yx = inpaint_nans(tmp0_yx, cfg.inpaint_method)
        tmp0_yx[np.isnan(tmp0_yx)] = 0
        nfert0 = np.nansum(
            test_croparea_x * 1e6 * out_cropfrac["garr_xv"][:, c2] * tmp0_yx.ravel(order="F")[np.asarray(gridlist["list2map"]) - 1] * 1e-3
        )
        nfert1 = np.sum(
            test_croparea_x * 1e6 * out_cropfrac["garr_xv"][:, c2] * tmp1_yx.ravel(order="F")[np.asarray(gridlist["list2map"]) - 1] * 1e-3
        )
        nfert_diff = nfert1 - nfert0
        nfert_diff_pct = nfert_diff / nfert0 * 100 if nfert0 != 0 else 0
        print(f"    Nfert increased {nfert_diff_pct:0.2f}% ({nfert_diff:0.1g} t/yr) in {interp_test_period_str}")
        out_nfert["garr_xv"][:, c2] = tmp1_yx.ravel(order="F")[np.asarray(gridlist["list2map"]) - 1]
        if np.isnan(out_nfert["garr_xv"][:, c2]).any():
            raise RuntimeError("NaN remaining in out_nfert.garr_xv(:,c2)!")

    if (out_cropfrac["garr_xv"] < 0).any():
        print("Warning: Setting negative members of out_cropfrac.garr_xv to zero.")
        out_cropfrac["garr_xv"][out_cropfrac["garr_xv"] < 0] = 0
    if (out_nfert["garr_xv"] < 0).any():
        print("Warning: Setting negative members of out_nfert.garr_xv to zero.")
        out_nfert["garr_xv"][out_nfert["garr_xv"] < 0] = 0

    tmp = np.sum(out_cropfrac["garr_xv"], axis=1)
    if np.any(np.abs(tmp - 1) > 1e-6):
        print("Warning: Normalizing cropfrac sums to 1")
        out_cropfrac["garr_xv"] = out_cropfrac["garr_xv"] / np.tile(tmp.reshape(-1, 1), (1, ncrops_out))

    # Add Miscanthus(i)
    out_cropfrac["garr_xv"] = np.column_stack([out_cropfrac["garr_xv"], np.zeros((n_cells, 2))])
    out_nfert["garr_xv"] = np.column_stack([out_nfert["garr_xv"], np.zeros((n_cells, 2))])
    out_cropfrac["varNames"] = list(out_cropfrac["varNames"]) + ["Miscanthus", "Miscanthusi"]
    out_nfert["varNames"] = list(out_nfert["varNames"]) + ["Miscanthus", "Miscanthusi"]

    # Headers and filenames
    out_lu_header = ["Lon", "Lat", "Year"] + list(out_lu["varNames"])
    out_cropfrac_header = ["Lon", "Lat"] + list(out_cropfrac["varNames"])
    out_nfert_header = ["Lon", "Lat"] + list(out_nfert["varNames"])

    out_file_lu = os.path.join(out_dir, f"LU.remapv{cfg.remap_ver}.txt")
    out_file_cropfrac = os.path.join(out_dir, f"cropfracs.remapv{cfg.remap_ver}.txt")
    out_file_nfert = os.path.join(out_dir, f"nfert.remapv{cfg.remap_ver}.txt")

    delete_existing_outputs(out_file_lu)
    delete_existing_outputs(out_file_cropfrac)
    delete_existing_outputs(out_file_nfert)

    # Save gridlist (random order)
    print("Saving gridlist...")
    lons4map = np.arange(-179.75, 180, 0.5)
    lats4map = np.arange(-89.75, 90, 0.5)
    lons_map = np.tile(lons4map, (len(lats4map), 1))
    lats_map = np.tile(lats4map.reshape(-1, 1), (1, len(lons4map)))
    lons_4gl = lons_map.ravel(order="F")[gridlist["mask_YX"].ravel(order="F")]
    lats_4gl = lats_map.ravel(order="F")[gridlist["mask_YX"].ravel(order="F")]
    n_cells_4gl = len(lons_4gl)
    rng = np.random.default_rng(20210106)
    rdmsam = rng.permutation(n_cells_4gl)
    out_file_gridlist = os.path.join(out_dir, f"gridlist.remapv{cfg.remap_ver}.txt")
    with open(out_file_gridlist, "w") as f:
        for lon, lat in zip(lons_4gl[rdmsam], lats_4gl[rdmsam]):
            f.write(f"{lon:4.2f} {lat:4.2f}\n")

    print("Saving LU...")
    save_table(
        out_lu_header,
        out_lu,
        out_file_lu,
        outPrec=cfg.out_prec,
        outWidth=cfg.out_width,
        delimiter=cfg.delimiter,
        overwrite=cfg.overwrite,
        fancy=cfg.fancy,
        progress_step_pct=1,
        gzip_output=cfg.do_gzip,
    )

    print("Saving cropfracs...")
    save_table(
        out_cropfrac_header,
        out_cropfrac,
        out_file_cropfrac,
        outPrec=cfg.out_prec,
        outWidth=cfg.out_width,
        delimiter=cfg.delimiter,
        overwrite=cfg.overwrite,
        fancy=cfg.fancy,
        progress_step_pct=1,
        gzip_output=cfg.do_gzip,
    )

    print("Saving nfert...")
    save_table(
        out_nfert_header,
        out_nfert,
        out_file_nfert,
        outPrec=cfg.out_prec,
        outWidth=cfg.out_width,
        delimiter=cfg.delimiter,
        overwrite=cfg.overwrite,
        fancy=cfg.fancy,
        progress_step_pct=1,
        gzip_output=cfg.do_gzip,
    )
    _restore_streams()


def _max_sum_diff(garr_xvy: np.ndarray, target: float) -> float:
    return float(np.max(np.abs(target - np.sum(garr_xvy, axis=1))))


class _TeeWriter:
    def __init__(self, primary: TextIO, secondary: TextIO) -> None:
        self.primary = primary
        self.secondary = secondary

    def write(self, data: str) -> int:
        self.primary.write(data)
        self.secondary.write(data)
        return len(data)

    def flush(self) -> None:
        self.primary.flush()
        self.secondary.flush()
