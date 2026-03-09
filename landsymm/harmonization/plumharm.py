"""Main PLUM harmonization driver (partial).

Pseudocode:
1) Load PLUMharm options + defaults.
2) Resolve plumDirs/harmDirs and validate.
3) Import baseline/reference data.
4) For each scenario and year:
   - Read PLUM inputs (LCF/Crop/Fert/Irrig).
   - Harmonize areas and management.
   - Enforce conservation checks.
   - Write yearly outputs.
"""

# TODO:
# - implement scenario/year loops with conservation checks
# - preserve all output naming conventions
from __future__ import annotations

import os
import time
from typing import List, Sequence

import h5py
import numpy as np
from scipy.io import savemat

from landsymm.common.aggregation import aggregate
from landsymm.common.lpjg_io import maps2table, save_table

from .plumharm_area import fix_tiny_negs, harmonize_area
from .plumharm_checks import check_bad_vals, check_mgmt_conservation
from .plumharm_debug import debug_out, debug_out_deltas
from .plumharm_dist import distribute_area_deltas, distribute_mgmt
from .plumharm_import_ref import import_ref_data
from .plumharm_mgmt import harmonize_mgmt
from .plumharm_options import PlumHarmConfig
from .plumharm_process_plum import process_plum_inputs


class PlumHarmPipeline:
    """Class-based PLUM harmonization pipeline (skeleton)."""

    def __init__(self, cfg: PlumHarmConfig) -> None:
        self.cfg = cfg

    def run(self) -> None:
        """Run the PLUM harmonization workflow (TODO)."""
        run_plumharm(self.cfg)


def run_plumharm(cfg: PlumHarmConfig) -> None:
    """Run PLUM harmonization (partial parity)."""
    if cfg.this_dir and not os.path.isdir(cfg.this_dir):
        raise RuntimeError(f"this_dir not found: {cfg.this_dir}")

    plum_dirs = _normalize_dirs(cfg.plum_dirs, cfg.this_dir)
    harm_dirs = cfg.harm_dirs
    if harm_dirs is None:
        harm_dirs = _get_harm_dirs(plum_dirs, cfg.fruitveg_sugar_2oil, cfg.combine_crops)
    harm_dirs = _normalize_dirs(harm_dirs, cfg.this_dir)

    if len(plum_dirs) != len(harm_dirs):
        raise RuntimeError("Numbers of plum_dirs and harm_dirs don't match")
    plum_dirs = _check_dirs(plum_dirs, "r")
    harm_dirs = _check_dirs(harm_dirs, "rw")
    for p, h in zip(plum_dirs, harm_dirs):
        if os.path.abspath(p) == os.path.abspath(h):
            raise RuntimeError(f"plumDir and harmDir must not be the same: {p}")

    if cfg.year1 <= cfg.base_year:
        raise RuntimeError("year1 must be greater than base_year")

    ref = import_ref_data(cfg, do_harm=True)
    land_area_yx = ref["landArea_YX"]
    land_area_2deg_yx = ref["landArea_2deg_YX"]
    base = ref["base"]
    base_2deg = ref["base_2deg"]
    base_2deg_bare_yx = ref["base_2deg_bare_YX"]
    lu_names = ref["lu_names"]
    lpjg_crops = ref["lpjg_crops"]
    plum_to_lpjg = ref["plum_to_lpjg"]
    not_bare = ref["notBare"]
    is_agri = ref["isAgri"]
    is_crop = ref["isCrop"]
    lu_names_agri = [v for v, ok in zip(lu_names, is_agri) if ok]
    is_agri_is_past = np.array([v == "PASTURE" for v in lu_names_agri])
    list2map = ref["list2map"]
    list2map_2deg = ref["list2map_2deg"]
    res_frac_terr_yx = ref["resFrac_terr_YX"]
    res_frac_prot_yx = ref["resFrac_prot_YX"]
    base_vegd_yx = ref["base_vegd_YX"]

    years = list(range(cfg.year1, cfg.yearN + 1))
    do_debug = cfg.debugIJ_2deg is not None and (
        cfg.debug_areas or cfg.debug_nfert or cfg.debug_irrig
    )
    save_halfdeg_any = cfg.save_halfdeg_mat or cfg.save_halfdeg_txt
    save_2deg_any = cfg.save_2deg_mat or cfg.save_2deg_txt
    save_any = save_halfdeg_any or save_2deg_any

    if cfg.combine_crops:
        print("Warning: Combining all crops into one! Will skip harmonization of nfert and irrig.")

    for plum_dir, harm_dir in zip(plum_dirs, harm_dirs):
        print(f"\n*\n*\n*\nPLUM output directory (plumDir): {plum_dir}")
        print(f"Harm. output directory (harmDir): {harm_dir}\n*\n*\n*")
        plum_config_file = os.path.join(plum_dir, "config.properties")
        if not os.path.isfile(plum_config_file):
            raise RuntimeError(f"PLUM config file not found: {plum_config_file}")
        min_natural_rate = _read_min_natural_rate(plum_config_file)
        res_frac = 1 - (1 - res_frac_terr_yx) * (1 - res_frac_prot_yx) * (1 - min_natural_rate)
        res_area_yx = res_frac * base_vegd_yx
        res_area_2deg_yx = aggregate(res_area_yx, 0.5, 2)
        if do_debug:
            db2i = cfg.debugIJ_2deg[0] - 1
            db2j = cfg.debugIJ_2deg[1] - 1
            lons_map_2deg = ref["lons_map_2deg"]
            lats_map_2deg = ref["lats_map_2deg"]
            print(
                "Center Lon {:.2f}, Lat {:.2f}, land area {:.4g}".format(
                    1 + lons_map_2deg[db2i, db2j],
                    1 + lats_map_2deg[db2i, db2j],
                    np.round(land_area_2deg_yx[db2i, db2j], 4),
                )
            )
            print("")
            if cfg.year1 == cfg.base_year + 1:
                if cfg.debug_areas:
                    debug_out(
                        "LU baseline",
                        "areas",
                        base_2deg["maps_YXv"],
                        land_area_2deg_yx,
                        cfg.debugIJ_2deg,
                        lu_names,
                    )
                if cfg.debug_nfert and not cfg.combine_crops:
                    debug_out(
                        "LU baseline",
                        "Nfert",
                        ref["base_nfert_tot_2deg"]["maps_YXv"],
                        base_2deg["maps_YXv"][:, :, is_crop],
                        cfg.debugIJ_2deg,
                        lpjg_crops,
                    )
                if cfg.debug_irrig and not cfg.combine_crops:
                    debug_out(
                        "LU baseline",
                        "irrig",
                        ref["base_irrig_tot_2deg"]["maps_YXv"],
                        base_2deg["maps_YXv"][:, :, is_crop],
                        cfg.debugIJ_2deg,
                        lpjg_crops,
                    )
        out_y0 = base
        out_y0_2deg = base_2deg
        out_y0_vegd_yx = np.sum(base["maps_YXv"][:, :, not_bare], axis=2)
        out_y0_2deg_vegd_yx = np.sum(base_2deg["maps_YXv"][:, :, not_bare], axis=2)
        out_y0_agri_yxv = base["maps_YXv"][:, :, is_agri]
        out_y0_2deg_agri_yxv = base_2deg["maps_YXv"][:, :, is_agri]
        out_y0_2deg_ntrl_yx = base_2deg["maps_YXv"][:, :, lu_names.index("NATURAL")]

        out_y0_nfert = ref["base_nfert"]["maps_YXv"] if not cfg.combine_crops else None
        out_y0_irrig = ref["base_irrig"]["maps_YXv"] if not cfg.combine_crops else None
        out_y0_2deg_nfert = (
            ref["base_2deg_nfert"]["maps_YXv"] if not cfg.combine_crops else None
        )
        out_y0_2deg_irrig = (
            ref["base_2deg_irrig"]["maps_YXv"] if not cfg.combine_crops else None
        )
        max_nfert_y0 = None
        max_orig_nfert_y0 = None

        latest_nfert = (
            -1 * np.ones((*land_area_2deg_yx.shape, len(lpjg_crops)))
            if cfg.use_latest_plum_mgmt and not cfg.combine_crops
            else None
        )
        latest_irrig = (
            -1 * np.ones((*land_area_2deg_yx.shape, len(lpjg_crops)))
            if cfg.use_latest_plum_mgmt and not cfg.combine_crops
            else None
        )

        in_y0 = None
        in_y0_2deg = None
        in_y0_nfert = None
        in_y0_irrig = None
        in_y0_2deg_nfert = None
        in_y0_2deg_irrig = None

        stop_years = False
        is_first_year_in_loop = True
        n_years = len(years)
        ts_in: dict = {}
        ts_out: dict = {}
        for this_year in years:
            print(str(this_year))
            tic = time.perf_counter()
            if (
                in_y0 is None
                and is_first_year_in_loop
                and this_year - 1 != cfg.base_year
            ):
                post_path = _find_post_restart_file(
                    harm_dir, this_year - 1, cfg.base_year
                )
                if post_path is not None:
                    if do_debug:
                        print(f"*y0* from {post_path}")
                    post = _load_post_mat(post_path)
                    in_y0 = post.get("in_y0")
                    in_y0_2deg = post.get("in_y0_2deg")
                    in_y0_agri = post.get("in_y0_agri_YXv")
                    in_y0_2deg_agri = post.get("in_y0_2deg_agri_YXv")
                    in_y0_nfert = post.get("in_y0_nfert")
                    in_y0_irrig = post.get("in_y0_irrig")
                    in_y0_2deg_nfert = post.get("in_y0_2deg_nfert")
                    in_y0_2deg_irrig = post.get("in_y0_2deg_irrig")
                    out_y0 = post.get("out_y0")
                    out_y0_2deg = post.get("out_y0_2deg")
                    out_y0_agri_yxv = post.get("out_y0_agri_YXv")
                    out_y0_2deg_agri_yxv = post.get("out_y0_2deg_agri_YXv")
                    out_y0_vegd_yx = post.get("out_y0_vegd_YX")
                    out_y0_2deg_vegd_yx = post.get("out_y0_2deg_vegd_YX")
                    out_y0_2deg_ntrl_yx = post.get("out_y0_2deg_ntrl_YX")
                    out_y0_nfert = post.get("out_y0_nfert_YXv")
                    out_y0_irrig = post.get("out_y0_irrig_YXv")
                    out_y0_2deg_nfert = _get_post_maps(
                        post, "out_y0_2deg_nfert", "out_y0_2deg_nfert_YXv"
                    )
                    out_y0_2deg_irrig = _get_post_maps(
                        post, "out_y0_2deg_irrig", "out_y0_2deg_irrig_YXv"
                    )
                    latest_nfert = post.get("latestPLUMin_2deg_nfert_YXv", latest_nfert)
                    latest_irrig = post.get("latestPLUMin_2deg_irrig_YXv", latest_irrig)
                    max_harm_nfert_y0 = post.get("max_harm_nfert_y0", None)
                    max_nfert_y0 = post.get("max_nfert_y0", max_nfert_y0)
                    max_orig_nfert_y0 = post.get("max_orig_nfert_y0", max_orig_nfert_y0)
                    bare_frac_y0 = post.get("bareFrac_y0_YX", None)
                    if in_y0 is not None and bare_frac_y0 is None:
                        bare_frac_y0 = (
                            in_y0["maps_YXv"][:, :, lu_names.index("BARREN")] / land_area_yx
                        )
                    if in_y0 is not None and in_y0_agri is None:
                        in_y0_agri = _align_maps_by_names(
                            in_y0["maps_YXv"], in_y0["varNames"], lu_names_agri
                        )
                    if in_y0_2deg is not None and in_y0_2deg_agri is None:
                        in_y0_2deg_agri = _align_maps_by_names(
                            in_y0_2deg["maps_YXv"], in_y0_2deg["varNames"], lu_names_agri
                        )
                    if out_y0 is not None and out_y0_agri_yxv is None:
                        out_y0_agri_yxv = out_y0["maps_YXv"][:, :, is_agri]
                    if out_y0_2deg is not None and out_y0_2deg_agri_yxv is None:
                        out_y0_2deg_agri_yxv = out_y0_2deg["maps_YXv"][:, :, is_agri]
                    if out_y0 is not None and out_y0_vegd_yx is None:
                        out_y0_vegd_yx = np.sum(out_y0["maps_YXv"][:, :, not_bare], axis=2)
                    if out_y0_2deg is not None and out_y0_2deg_vegd_yx is None:
                        out_y0_2deg_vegd_yx = np.sum(
                            out_y0_2deg["maps_YXv"][:, :, not_bare], axis=2
                        )
                    if out_y0_2deg is not None and out_y0_2deg_ntrl_yx is None:
                        out_y0_2deg_ntrl_yx = out_y0_2deg["maps_YXv"][
                            :, :, lu_names.index("NATURAL")
                        ]
            if in_y0 is None:
                file_in = os.path.join(plum_dir, str(this_year - 1), "LandCoverFract.txt")
                if do_debug and this_year - 1 == cfg.base_year:
                    print("out_y0 from baseline LU")
                if do_debug:
                    print(f"in_y0 from {file_in}")
                (
                    in_y0,
                    in_y0_nfert,
                    in_y0_irrig,
                    in_y0_2deg,
                    in_y0_2deg_nfert,
                    in_y0_2deg_irrig,
                    latest_nfert,
                    latest_irrig,
                    max_orig_nfert_y0,
                ) = process_plum_inputs(
                    file_in,
                    land_area_yx,
                    land_area_2deg_yx,
                    lu_names,
                    None,
                    latest_nfert,
                    latest_irrig,
                    plum_to_lpjg,
                    lpjg_crops,
                    cfg.norm2extra,
                    cfg.inpaint_method,
                    cfg.fruitveg_sugar_2oil,
                    cfg.allow_unveg,
                )
                bare_frac_y0 = (
                    in_y0["maps_YXv"][:, :, lu_names.index("BARREN")] / land_area_yx
                )
                check_bad_vals(
                    in_y0["maps_YXv"],
                    None,
                    None,
                    land_area_yx,
                    lu_names,
                    "in_y0",
                    cfg.out_prec,
                )

            file_in = os.path.join(plum_dir, str(this_year), "LandCoverFract.txt")
            if do_debug:
                print(f"in_y1 from {file_in}")
            (
                in_y1,
                in_y1_nfert,
                in_y1_irrig,
                in_y1_2deg,
                in_y1_2deg_nfert,
                in_y1_2deg_irrig,
                latest_nfert,
                latest_irrig,
                max_orig_nfert_y1,
            ) = process_plum_inputs(
                file_in,
                land_area_yx,
                land_area_2deg_yx,
                lu_names,
                bare_frac_y0,
                latest_nfert,
                latest_irrig,
                plum_to_lpjg,
                lpjg_crops,
                cfg.norm2extra,
                cfg.inpaint_method,
                cfg.fruitveg_sugar_2oil,
                cfg.allow_unveg,
            )
            check_bad_vals(
                in_y1["maps_YXv"],
                None,
                None,
                land_area_yx,
                lu_names,
                "in_y1",
                cfg.out_prec,
            )

            in_y0_agri = _align_maps_by_names(
                in_y0["maps_YXv"], in_y0["varNames"], lu_names_agri
            )
            in_y1_agri = _align_maps_by_names(
                in_y1["maps_YXv"], in_y1["varNames"], lu_names_agri
            )
            in_y0_2deg_agri = _align_maps_by_names(
                in_y0_2deg["maps_YXv"], in_y0_2deg["varNames"], lu_names_agri
            )
            in_y1_2deg_agri = _align_maps_by_names(
                in_y1_2deg["maps_YXv"], in_y1_2deg["varNames"], lu_names_agri
            )

            if this_year - 1 == cfg.base_year and not cfg.combine_crops:
                print(f"    PLUM baseline ({cfg.base_year}) discrepancies relative to baseline LU:")
                for i, this_var in enumerate(in_y0["varNames"]):
                    tmp_in = np.nansum(in_y0["maps_YXv"][:, :, i])
                    tmp_out = np.nansum(base["maps_YXv"][:, :, i])
                    _check_discrepancy(tmp_in, tmp_out, this_var, "area")
                    if not is_crop[i]:
                        continue
                    i_crop = lpjg_crops.index(this_var)
                    tmp_in = np.nansum(
                        in_y0_nfert["maps_YXv"][:, :, i_crop] * in_y0["maps_YXv"][:, :, i]
                    )
                    tmp_out = np.nansum(
                        out_y0_nfert[:, :, i_crop] * base["maps_YXv"][:, :, i]
                    )
                    _check_discrepancy(tmp_in, tmp_out, this_var, "fert")
                    tmp_in = np.nansum(
                        in_y0_irrig["maps_YXv"][:, :, i_crop] * in_y0["maps_YXv"][:, :, i]
                    )
                    tmp_out = np.nansum(
                        out_y0_irrig[:, :, i_crop] * base["maps_YXv"][:, :, i]
                    )
                    _check_discrepancy(tmp_in, tmp_out, this_var, "irrig")
            if cfg.debug_areas and not cfg.combine_crops:
                _debug_global_areas(
                    in_y0_2deg["maps_YXv"],
                    in_y1_2deg["maps_YXv"],
                    "Initial import",
                    "in_2deg",
                    "in_2deg",
                    lu_names,
                    is_crop,
                    is_crop,
                    cfg.dbCrop,
                    this_year,
                )
            if cfg.debug_nfert and not cfg.combine_crops:
                _debug_global_mgmts(
                    in_y0_2deg_nfert["maps_YXv"],
                    in_y1_2deg_nfert["maps_YXv"],
                    np.zeros_like(in_y1_2deg_nfert["maps_YXv"]),
                    in_y0_2deg_agri[:, :, ~is_agri_is_past],
                    in_y1_2deg_agri[:, :, ~is_agri_is_past],
                    "Mt",
                    1e-6 * 1e-3,
                    "Initial import",
                    "in_2deg",
                    "in_2deg",
                    "nfert",
                    lpjg_crops,
                    cfg.dbCrop,
                    this_year,
                )
            if cfg.debug_irrig and not cfg.combine_crops:
                _debug_global_mgmts(
                    in_y0_2deg_irrig["maps_YXv"],
                    in_y1_2deg_irrig["maps_YXv"],
                    np.zeros_like(in_y1_2deg_irrig["maps_YXv"]),
                    in_y0_2deg_agri[:, :, ~is_agri_is_past],
                    in_y1_2deg_agri[:, :, ~is_agri_is_past],
                    "arb",
                    1,
                    "Initial import",
                    "in_2deg",
                    "in_2deg",
                    "irrig",
                    lpjg_crops,
                    cfg.dbCrop,
                    this_year,
                )

            set_area_to_zero, updated_in_y0 = _get_set_area_to_zero(
                lu_names,
                in_y0["maps_YXv"],
                in_y1["maps_YXv"],
                out_y0["maps_YXv"] if isinstance(out_y0, dict) else base["maps_YXv"],
                this_year,
                cfg.base_year,
                cfg.conserv_tol_area,
            )
            if updated_in_y0:
                in_y0_2deg["maps_YXv"] = aggregate(in_y0["maps_YXv"], 0.5, 2)
                in_y0_agri = _align_maps_by_names(
                    in_y0["maps_YXv"], in_y0["varNames"], lu_names_agri
                )
                in_y0_2deg_agri = _align_maps_by_names(
                    in_y0_2deg["maps_YXv"], in_y0_2deg["varNames"], lu_names_agri
                )
            check_bad_vals(
                in_y0["maps_YXv"],
                None,
                None,
                land_area_yx,
                lu_names,
                "in_y0.2",
                cfg.out_prec,
            )
            check_bad_vals(
                in_y1["maps_YXv"],
                None,
                None,
                land_area_yx,
                lu_names,
                "in_y1.2",
                cfg.out_prec,
            )
            if cfg.debug_areas and cfg.debugIJ_2deg is not None:
                debug_out_deltas(
                    "iny0_to_iny1",
                    "areas",
                    in_y0_2deg["maps_YXv"],
                    in_y1_2deg["maps_YXv"],
                    cfg.debugIJ_2deg,
                    lu_names,
                )

            # MATLAB order parity: refresh current y0 NATURAL from out_y0_2deg
            # maps each year before area harmonization.
            out_y0_2deg_ntrl_yx = out_y0_2deg["maps_YXv"][:, :, lu_names.index("NATURAL")]

            try:
                (
                    out_y1_2deg_agri,
                    out_y1_2deg_ntrl,
                    out_y1_2deg_bare,
                ) = harmonize_area(
                    in_y0_2deg_agri,
                    in_y1_2deg_agri,
                    out_y0_2deg_agri_yxv,
                    out_y0_2deg_vegd_yx,
                    out_y0_2deg["maps_YXv"][:, :, lu_names.index("NATURAL")],
                    base_2deg_bare_yx,
                    land_area_2deg_yx,
                    res_area_2deg_yx,
                    lu_names,
                    lu_names_agri,
                    is_agri,
                    set_area_to_zero,
                    cfg.fix_tiny_negs_tol_m2,
                    cfg.conserv_tol_pct,
                    cfg.conserv_tol_area,
                    cfg.debugIJ_2deg,
                    cfg.dbCrop,
                )
            except RuntimeError as exc:
                msg = str(exc)
                if "Possible infinite loop" in msg:
                    print("INFINITE LOOP IN RINGREDIST AREA; HALTING.")
                    stop_years = True
                else:
                    raise
            if stop_years:
                break

            out_y1_agri = distribute_area_deltas(
                land_area_yx,
                land_area_2deg_yx,
                out_y0_2deg_agri_yxv,
                out_y1_2deg_agri,
                out_y0_agri_yxv,
                out_y0_vegd_yx,
                cfg.fix_tiny_negs_tol_m2,
                cfg.conserv_tol_pct,
                cfg.conserv_tol_area,
                lu_names_agri,
                None,
            )
            if cfg.debugIJ_2deg is not None and cfg.dbCrop:
                try:
                    db_crop_idx = (
                        int(cfg.dbCrop)
                        if str(cfg.dbCrop).isdigit()
                        else lu_names_agri.index(cfg.dbCrop)
                    )
                except (ValueError, IndexError):
                    db_crop_idx = None
                if db_crop_idx is not None:
                    i2 = cfg.debugIJ_2deg[0] - 1
                    j2 = cfg.debugIJ_2deg[1] - 1
                    y0 = 4 * i2
                    x0 = 4 * j2
                    half_sum = float(
                        np.sum(out_y1_agri[y0 : y0 + 4, x0 : x0 + 4, db_crop_idx])
                    )
                    two_deg = float(out_y1_2deg_agri[i2, j2, db_crop_idx])
                    print(
                        "        DEBUG_AREA_HALF_SUM:",
                        f"ij=({cfg.debugIJ_2deg[0]},{cfg.debugIJ_2deg[1]})",
                        f"crop={lu_names_agri[db_crop_idx]}",
                        f"half_sum={half_sum:.12e}",
                        f"two_deg={two_deg:.12e}",
                    )
            if cfg.debug_areas:
                _debug_global_areas(
                    out_y0_2deg["maps_YXv"],
                    np.concatenate(
                        [
                            out_y1_2deg_agri,
                            out_y1_2deg_ntrl[:, :, None],
                            out_y1_2deg_bare[:, :, None],
                        ],
                        axis=2,
                    ),
                    "After ringRedist",
                    "out_2deg",
                    "out_2deg",
                    lu_names,
                    is_crop,
                    is_agri,
                    cfg.dbCrop,
                    this_year,
                )
                if cfg.debugIJ_2deg is not None:
                    debug_out_deltas(
                        "outy0_to_outy1",
                        "areas",
                        out_y0_2deg["maps_YXv"],
                        np.concatenate(
                            [
                                out_y1_2deg_agri,
                                out_y1_2deg_ntrl[:, :, None],
                                out_y1_2deg_bare[:, :, None],
                            ],
                            axis=2,
                        ),
                        cfg.debugIJ_2deg,
                        lu_names,
                    )
            _check_preserved_global_deltas(
                out_y0_2deg_agri_yxv,
                out_y1_2deg_agri,
                out_y0_agri_yxv,
                out_y1_agri,
                cfg.conserv_tol_pct,
                lu_names_agri,
                "2-deg to half-deg",
                "agri",
            )
            out_y1_ntrl = out_y0_vegd_yx - np.sum(out_y1_agri, axis=2)
            out_y1_vegd = np.sum(out_y1_agri, axis=2) + out_y1_ntrl
            out_y1_bare = land_area_yx - out_y1_vegd
            if cfg.debug_areas:
                _debug_global_areas(
                    out_y0_2deg["maps_YXv"],
                    np.concatenate(
                        [
                            out_y1_agri,
                            out_y1_ntrl[:, :, None],
                            out_y1_bare[:, :, None],
                        ],
                        axis=2,
                    ),
                    "Now at half-degree but before negative value fix",
                    "out",
                    "out",
                    lu_names,
                    is_crop,
                    is_agri,
                    cfg.dbCrop,
                    this_year,
                )
            tmp = np.concatenate(
                [out_y1_agri, out_y1_ntrl[:, :, None], out_y1_bare[:, :, None]], axis=2
            )
            tmp = fix_tiny_negs(
                tmp,
                np.repeat(land_area_yx[:, :, None], tmp.shape[2], axis=2),
                lu_names,
                cfg.out_prec,
                cfg.fix_tiny_negs_tol_m2,
                cfg.conserv_tol_area,
                None,
            )
            out_y1_agri = tmp[:, :, : -2]
            out_y1_ntrl = tmp[:, :, -2]
            out_y1_bare = tmp[:, :, -1]
            if cfg.debug_areas:
                _debug_global_areas(
                    out_y0_2deg["maps_YXv"],
                    np.concatenate(
                        [
                            out_y1_agri,
                            out_y1_ntrl[:, :, None],
                            out_y1_bare[:, :, None],
                        ],
                        axis=2,
                    ),
                    "At half-degree after negative value fix",
                    "out",
                    "out",
                    lu_names,
                    is_crop,
                    is_agri,
                    cfg.dbCrop,
                    this_year,
                )
            check_bad_vals(
                np.concatenate(
                    [out_y1_agri, out_y1_ntrl[:, :, None], out_y1_bare[:, :, None]], axis=2
                ),
                None,
                None,
                land_area_yx,
                lu_names,
                "out_y1",
                cfg.out_prec,
            )

            out_y1_nfert = None
            out_y1_irrig = None
            out_y1_2deg_nfert = None
            out_y1_2deg_irrig = None

            if not cfg.combine_crops:
                max_harm_nfert_y0 = np.nanmax(out_y0_2deg_nfert, axis=(0, 1))
                max_stack = [max_orig_nfert_y0, max_orig_nfert_y1, max_harm_nfert_y0]
                if max_nfert_y0 is not None:
                    max_stack.append(max_nfert_y0)
                max_nfert_y1 = np.max(np.vstack(max_stack), axis=0)
                max_nfert_y1_yxv = np.broadcast_to(
                    max_nfert_y1, in_y0_2deg_nfert["maps_YXv"].shape
                )
                too_much_nfert = in_y0_2deg_nfert["maps_YXv"] > max_nfert_y1_yxv
                in_y0_2deg_nfert["maps_YXv"][in_y0_2deg_nfert["maps_YXv"] < 0] = 0
                in_y1_2deg_nfert["maps_YXv"][in_y1_2deg_nfert["maps_YXv"] < 0] = 0
                in_y0_2deg_nfert["maps_YXv"][too_much_nfert] = max_nfert_y1_yxv[too_much_nfert]
                in_y1_2deg_nfert["maps_YXv"][too_much_nfert] = max_nfert_y1_yxv[too_much_nfert]
                in_y0_2deg_irrig["maps_YXv"][in_y0_2deg_irrig["maps_YXv"] < 0] = 0
                in_y1_2deg_irrig["maps_YXv"][in_y1_2deg_irrig["maps_YXv"] < 0] = 0
                in_y0_2deg_irrig["maps_YXv"][in_y0_2deg_irrig["maps_YXv"] > 1] = 1
                in_y1_2deg_irrig["maps_YXv"][in_y1_2deg_irrig["maps_YXv"] > 1] = 1

                check_bad_vals(
                    None,
                    in_y0_2deg_nfert["maps_YXv"],
                    in_y0_2deg_irrig["maps_YXv"],
                    None,
                    lu_names,
                    "in_y0_2deg",
                    cfg.out_prec,
                )
                check_bad_vals(
                    None,
                    in_y1_2deg_nfert["maps_YXv"],
                    in_y1_2deg_irrig["maps_YXv"],
                    None,
                    lu_names,
                    "in_y1_2deg",
                    cfg.out_prec,
                )
                if cfg.debug_nfert:
                    _debug_global_mgmts(
                        in_y0_2deg_nfert["maps_YXv"],
                        in_y1_2deg_nfert["maps_YXv"],
                        np.zeros_like(in_y1_2deg_nfert["maps_YXv"]),
                        in_y0_2deg_agri[:, :, ~is_agri_is_past],
                        in_y1_2deg_agri[:, :, ~is_agri_is_past],
                        "Mt",
                        1e-6 * 1e-3,
                        "After limiting to [0, max]",
                        "in_2deg",
                        "in_2deg",
                        "nfert",
                        lpjg_crops,
                        cfg.dbCrop,
                        this_year,
                    )
                    if cfg.debugIJ_2deg is not None:
                        debug_out_deltas(
                            "iny0_to_iny1",
                            "Nfert",
                            in_y0_2deg_nfert["maps_YXv"] * in_y0_2deg_agri[:, :, ~is_agri_is_past],
                            in_y1_2deg_nfert["maps_YXv"] * in_y1_2deg_agri[:, :, ~is_agri_is_past],
                            cfg.debugIJ_2deg,
                            lpjg_crops,
                        )
                if cfg.debug_irrig:
                    _debug_global_mgmts(
                        in_y0_2deg_irrig["maps_YXv"],
                        in_y1_2deg_irrig["maps_YXv"],
                        np.zeros_like(in_y1_2deg_irrig["maps_YXv"]),
                        in_y0_2deg_agri[:, :, ~is_agri_is_past],
                        in_y1_2deg_agri[:, :, ~is_agri_is_past],
                        "arb",
                        1,
                        "After limiting to [0, max]",
                        "in_2deg",
                        "in_2deg",
                        "irrig",
                        lpjg_crops,
                        cfg.dbCrop,
                        this_year,
                    )
                    if cfg.debugIJ_2deg is not None:
                        debug_out_deltas(
                            "iny0_to_iny1",
                            "irrig",
                            in_y0_2deg_irrig["maps_YXv"] * in_y0_2deg_agri[:, :, ~is_agri_is_past],
                            in_y1_2deg_irrig["maps_YXv"] * in_y1_2deg_agri[:, :, ~is_agri_is_past],
                            cfg.debugIJ_2deg,
                            lpjg_crops,
                        )
                try:
                    (
                        out_y1_2deg_nfert,
                        out_y1_2deg_irrig,
                        _unm2_nfert,
                        _unm2_irrig,
                        not_enough_nfert,
                        not_enough_irrig,
                        mid_y1_2deg_nfert,
                        mid_y1_2deg_irrig,
                        _unm_y1_2deg_nfert,
                        _unm_y1_2deg_irrig,
                    ) = harmonize_mgmt(
                        out_y0_2deg_nfert,
                        out_y0_2deg_irrig,
                        out_y0_2deg_agri_yxv[:, :, ~is_agri_is_past],
                        in_y0_2deg_nfert["maps_YXv"],
                        in_y1_2deg_nfert["maps_YXv"],
                        in_y0_2deg_irrig["maps_YXv"],
                        in_y1_2deg_irrig["maps_YXv"],
                        in_y0_2deg_agri[:, :, ~is_agri_is_past],
                        in_y1_2deg_agri[:, :, ~is_agri_is_past],
                        out_y1_2deg_agri[:, :, ~is_agri_is_past],
                        max_nfert_y1,
                        lpjg_crops,
                        cfg.conserv_tol_pct,
                        cfg.debugIJ_2deg,
                        cfg.dbCrop,
                        cfg.debug_nfert,
                        cfg.debug_irrig,
                        cfg.out_prec,
                    )
                except RuntimeError as exc:
                    msg = str(exc)
                    if "Possible infinite loop" in msg and "ringRedist nfert" in msg:
                        print("POSSIBLE INFINITE LOOP IN RINGREDIST NFERT; HALTING.")
                        stop_years = True
                    elif "Possible infinite loop" in msg and "ringRedist irrig" in msg:
                        print("POSSIBLE INFINITE LOOP IN RINGREDIST IRRIG; HALTING.")
                        stop_years = True
                    else:
                        raise
                if stop_years:
                    break
                if cfg.debug_nfert:
                    _debug_global_mgmts(
                        out_y0_2deg_nfert,
                        mid_y1_2deg_nfert,
                        _unm_y1_2deg_nfert,
                        out_y0_2deg_agri_yxv[:, :, ~is_agri_is_past],
                        out_y1_2deg_agri[:, :, ~is_agri_is_past],
                        "Mt",
                        1e-6 * 1e-3,
                        "After applying deltas",
                        "out_2deg",
                        "mid_2deg",
                        "nfert",
                        lpjg_crops,
                        cfg.dbCrop,
                        this_year,
                    )
                if cfg.debug_irrig:
                    _debug_global_mgmts(
                        out_y0_2deg_irrig,
                        mid_y1_2deg_irrig,
                        _unm_y1_2deg_irrig,
                        out_y0_2deg_agri_yxv[:, :, ~is_agri_is_past],
                        out_y1_2deg_agri[:, :, ~is_agri_is_past],
                        "arb",
                        1,
                        "After applying deltas",
                        "out_2deg",
                        "mid_2deg",
                        "irrig",
                        lpjg_crops,
                        cfg.dbCrop,
                        this_year,
                    )
                if cfg.debug_nfert:
                    _debug_global_mgmts(
                        out_y0_2deg_nfert,
                        out_y1_2deg_nfert,
                        np.zeros_like(out_y1_2deg_nfert),
                        out_y0_2deg_agri_yxv[:, :, ~is_agri_is_past],
                        out_y1_2deg_agri[:, :, ~is_agri_is_past],
                        "Mt",
                        1e-6 * 1e-3,
                        "After ringRedist",
                        "out_2deg",
                        "out_2deg",
                        "nfert",
                        lpjg_crops,
                        cfg.dbCrop,
                        this_year,
                    )
                if cfg.debug_irrig:
                    _debug_global_mgmts(
                        out_y0_2deg_irrig,
                        out_y1_2deg_irrig,
                        np.zeros_like(out_y1_2deg_irrig),
                        out_y0_2deg_agri_yxv[:, :, ~is_agri_is_past],
                        out_y1_2deg_agri[:, :, ~is_agri_is_past],
                        "arb",
                        1,
                        "After ringRedist",
                        "out_2deg",
                        "out_2deg",
                        "irrig",
                        lpjg_crops,
                        cfg.dbCrop,
                        this_year,
                    )
                if cfg.debug_nfert and cfg.debugIJ_2deg is not None and cfg.dbCrop is not None:
                    try:
                        db_idx = int(cfg.dbCrop)
                    except (TypeError, ValueError):
                        db_idx = (
                            list(lpjg_crops).index(cfg.dbCrop)
                            if cfg.dbCrop in list(lpjg_crops)
                            else None
                        )
                    if db_idx is not None:
                        i = cfg.debugIJ_2deg[0] - 1
                        j = cfg.debugIJ_2deg[1] - 1
                        if out_y1_2deg_nfert is not None and out_y1_2deg_agri is not None:
                            rate = float(out_y1_2deg_nfert[i, j, db_idx])
                            area = float(out_y1_2deg_agri[i, j, db_idx])
                            total = rate * area
                            print(
                                "        DEBUG_OUTPUT_NFERT_2DEG:",
                                f"ij=({cfg.debugIJ_2deg[0]},{cfg.debugIJ_2deg[1]})",
                                f"crop={lpjg_crops[db_idx]}",
                                f"rate={rate:.17e}",
                                f"area={area:.17e}",
                                f"total={total:.17e}",
                            )
                        if out_y1_2deg_irrig is not None and out_y1_2deg_agri is not None:
                            rate = float(out_y1_2deg_irrig[i, j, db_idx])
                            area = float(out_y1_2deg_agri[i, j, db_idx])
                            total = rate * area
                            print(
                                "        DEBUG_OUTPUT_IRRIG_2DEG:",
                                f"ij=({cfg.debugIJ_2deg[0]},{cfg.debugIJ_2deg[1]})",
                                f"crop={lpjg_crops[db_idx]}",
                                f"rate={rate:.17e}",
                                f"area={area:.17e}",
                                f"total={total:.17e}",
                            )
                check_mgmt_conservation(
                    out_y0_2deg_nfert,
                    out_y0_2deg_agri_yxv[:, :, ~is_agri_is_past],
                    out_y1_2deg_nfert,
                    out_y1_2deg_agri[:, :, ~is_agri_is_past],
                    in_y0_2deg_nfert["maps_YXv"],
                    in_y0_2deg_agri[:, :, ~is_agri_is_past],
                    in_y1_2deg_nfert["maps_YXv"],
                    in_y1_2deg_agri[:, :, ~is_agri_is_past],
                    _unm2_nfert,
                    lpjg_crops,
                    cfg.conserv_tol_pct,
                    not_enough_nfert,
                    "3b nfert",
                    True,
                )
                check_mgmt_conservation(
                    out_y0_2deg_irrig,
                    out_y0_2deg_agri_yxv[:, :, ~is_agri_is_past],
                    out_y1_2deg_irrig,
                    out_y1_2deg_agri[:, :, ~is_agri_is_past],
                    in_y0_2deg_irrig["maps_YXv"],
                    in_y0_2deg_agri[:, :, ~is_agri_is_past],
                    in_y1_2deg_irrig["maps_YXv"],
                    in_y1_2deg_agri[:, :, ~is_agri_is_past],
                    _unm2_irrig,
                    lpjg_crops,
                    cfg.conserv_tol_pct,
                    not_enough_irrig,
                    "3b irrig",
                    True,
                )
                out_y1_nfert = distribute_mgmt(out_y1_2deg_nfert, 2, 0.5)
                out_y1_irrig = distribute_mgmt(out_y1_2deg_irrig, 2, 0.5)
                out_y1_nfert[out_y1_agri[:, :, ~is_agri_is_past] == 0] = 0
                out_y1_irrig[out_y1_agri[:, :, ~is_agri_is_past] == 0] = 0
                check_bad_vals(
                    None,
                    out_y1_nfert,
                    out_y1_irrig,
                    None,
                    lu_names,
                    "out_y1",
                    cfg.out_prec,
                )

                if cfg.debug_nfert:
                    _debug_global_mgmts(
                        out_y0_nfert,
                        out_y1_nfert,
                        np.zeros_like(out_y1_nfert),
                        out_y0_agri_yxv[:, :, ~is_agri_is_past],
                        out_y1_agri[:, :, ~is_agri_is_past],
                        "Mt",
                        1e-6 * 1e-3,
                        "Now at half-degree",
                        "out",
                        "out",
                        "nfert",
                        lpjg_crops,
                        cfg.dbCrop,
                        this_year,
                    )
                if cfg.debug_irrig:
                    _debug_global_mgmts(
                        out_y0_irrig,
                        out_y1_irrig,
                        np.zeros_like(out_y1_irrig),
                        out_y0_agri_yxv[:, :, ~is_agri_is_past],
                        out_y1_agri[:, :, ~is_agri_is_past],
                        "arb",
                        1,
                        "Now at half-degree",
                        "out",
                        "out",
                        "irrig",
                        lpjg_crops,
                        cfg.dbCrop,
                        this_year,
                    )
                check_mgmt_conservation(
                    out_y0_nfert,
                    out_y0_agri_yxv[:, :, ~is_agri_is_past],
                    out_y1_nfert,
                    out_y1_agri[:, :, ~is_agri_is_past],
                    in_y0_2deg_nfert["maps_YXv"],
                    in_y0_2deg_agri[:, :, ~is_agri_is_past],
                    in_y1_2deg_nfert["maps_YXv"],
                    in_y1_2deg_agri[:, :, ~is_agri_is_past],
                    np.zeros_like(out_y0_nfert),
                    lpjg_crops,
                    cfg.conserv_tol_pct,
                    not_enough_nfert,
                    "5 nfert",
                    False,
                )
                check_mgmt_conservation(
                    out_y0_irrig,
                    out_y0_agri_yxv[:, :, ~is_agri_is_past],
                    out_y1_irrig,
                    out_y1_agri[:, :, ~is_agri_is_past],
                    in_y0_2deg_irrig["maps_YXv"],
                    in_y0_2deg_agri[:, :, ~is_agri_is_past],
                    in_y1_2deg_irrig["maps_YXv"],
                    in_y1_2deg_agri[:, :, ~is_agri_is_past],
                    np.zeros_like(out_y0_irrig),
                    lpjg_crops,
                    cfg.conserv_tol_pct,
                    not_enough_irrig,
                    "5 irrig",
                    False,
                )
                if "ExtraCrop" in lpjg_crops:
                    ix_extra = lpjg_crops.index("ExtraCrop")
                    if np.max(out_y1_nfert[:, :, ix_extra]) > 0:
                        raise RuntimeError(
                            f"Some ExtraCrop Nfert >0! ({np.max(out_y1_nfert[:, :, ix_extra]):0.3e})"
                        )
                    if np.max(out_y1_irrig[:, :, ix_extra]) > 0:
                        raise RuntimeError(
                            f"Some ExtraCrop irrig >0! ({np.max(out_y1_irrig[:, :, ix_extra]):0.3e})"
                        )

            if cfg.save_halfdeg_mat or cfg.save_2deg_mat:
                harm_dir_year = os.path.join(harm_dir, str(this_year))
                os.makedirs(harm_dir_year, exist_ok=True)

            if save_any:
                print(f"  Done processing ({_toc_hms(time.perf_counter() - tic)}). Now writing.")

            if cfg.save_halfdeg_mat:
                out_lcf = {
                    "varNames": np.array(["PASTURE", "CROPLAND", "NATURAL", "BARREN"], dtype=object),
                    "maps_YXv": np.stack(
                        [
                            out_y1_agri[:, :, is_agri_is_past].squeeze(axis=2),
                            np.sum(out_y1_agri[:, :, ~is_agri_is_past], axis=2),
                            out_y1_ntrl,
                            out_y1_bare,
                        ],
                        axis=2,
                    )
                    / np.repeat(land_area_yx[:, :, None], 4, axis=2),
                }
                savemat(
                    os.path.join(harm_dir_year, f"LandCoverFract.base{cfg.base_year}.mat"),
                    {"out_y1": out_lcf},
                )
                if not cfg.combine_crops:
                    denom = np.sum(out_y1_agri[:, :, ~is_agri_is_past], axis=2)
                    out_cf = {
                        "varNames": np.array(list(lpjg_crops), dtype=object),
                        "maps_YXv": out_y1_agri[:, :, ~is_agri_is_past]
                        / np.repeat(denom[:, :, None], len(lpjg_crops), axis=2),
                    }
                    out_cf["maps_YXv"][np.isnan(out_cf["maps_YXv"])] = 0
                    out_cf["maps_YXv"][
                        np.repeat(denom[:, :, None] == 0, len(lpjg_crops), axis=2)
                    ] = 0
                    savemat(
                        os.path.join(harm_dir_year, f"CropFract.base{cfg.base_year}.mat"),
                        {"out_y1": out_cf},
                    )
                    out_fert = {
                        "varNames": np.array(list(lpjg_crops), dtype=object),
                        "maps_YXv": 1e4 * out_y1_nfert,
                    }
                    savemat(
                        os.path.join(harm_dir_year, f"Fert.base{cfg.base_year}.mat"),
                        {"out_y1": out_fert},
                    )
                    out_ir = {
                        "varNames": np.array(list(lpjg_crops), dtype=object),
                        "maps_YXv": out_y1_irrig,
                    }
                    savemat(
                        os.path.join(harm_dir_year, f"Irrig.base{cfg.base_year}.mat"),
                        {"out_y1": out_ir},
                    )
            if cfg.save_halfdeg_txt:
                out_lcf_txt = {
                    "varNames": ["PASTURE", "CROPLAND", "NATURAL", "BARREN"],
                    "maps_YXv": np.stack(
                        [
                            out_y1_agri[:, :, is_agri_is_past].squeeze(axis=2),
                            np.sum(out_y1_agri[:, :, ~is_agri_is_past], axis=2),
                            out_y1_ntrl,
                            out_y1_bare,
                        ],
                        axis=2,
                    )
                    / np.repeat(land_area_yx[:, :, None], 4, axis=2),
                }
                out_array, out_header = maps2table(out_lcf_txt, list2map)
                save_table(
                    out_header,
                    out_array,
                    os.path.join(harm_dir_year, f"LandCoverFract.base{cfg.base_year}.txt"),
                    outPrec=cfg.out_prec,
                    outWidth=cfg.out_width,
                    delimiter=cfg.delimiter,
                    overwrite=cfg.overwrite,
                    fancy=cfg.fancy,
                    progress_step_pct=20,
                    verbose=False,
                )
                if not cfg.combine_crops:
                    denom = np.sum(out_y1_agri[:, :, ~is_agri_is_past], axis=2)
                    out_cf_txt = {
                        "varNames": list(lpjg_crops),
                        "maps_YXv": out_y1_agri[:, :, ~is_agri_is_past]
                        / np.repeat(denom[:, :, None], len(lpjg_crops), axis=2),
                    }
                    out_cf_txt["maps_YXv"][np.isnan(out_cf_txt["maps_YXv"])] = 0
                    out_array, out_header = maps2table(out_cf_txt, list2map)
                    save_table(
                        out_header,
                        out_array,
                        os.path.join(harm_dir_year, f"CropFract.base{cfg.base_year}.txt"),
                        outPrec=cfg.out_prec,
                        outWidth=cfg.out_width,
                        delimiter=cfg.delimiter,
                        overwrite=cfg.overwrite,
                        fancy=cfg.fancy,
                        progress_step_pct=20,
                        verbose=False,
                    )
                    out_fert_txt = {"varNames": list(lpjg_crops), "maps_YXv": 1e4 * out_y1_nfert}
                    out_array, out_header = maps2table(out_fert_txt, list2map)
                    save_table(
                        out_header,
                        out_array,
                        os.path.join(harm_dir_year, f"Fert.base{cfg.base_year}.txt"),
                        outPrec=cfg.out_prec,
                        outWidth=cfg.out_width,
                        delimiter=cfg.delimiter,
                        overwrite=cfg.overwrite,
                        fancy=cfg.fancy,
                        progress_step_pct=20,
                        verbose=False,
                    )
                    out_ir_txt = {"varNames": list(lpjg_crops), "maps_YXv": out_y1_irrig}
                    out_array, out_header = maps2table(out_ir_txt, list2map)
                    save_table(
                        out_header,
                        out_array,
                        os.path.join(harm_dir_year, f"Irrig.base{cfg.base_year}.txt"),
                        outPrec=cfg.out_prec,
                        outWidth=cfg.out_width,
                        delimiter=cfg.delimiter,
                        overwrite=cfg.overwrite,
                        fancy=cfg.fancy,
                        progress_step_pct=20,
                        verbose=False,
                    )

            if cfg.save_2deg_mat:
                out_lcf_2 = {
                    "varNames": np.array(["PASTURE", "CROPLAND", "NATURAL", "BARREN"], dtype=object),
                    "maps_YXv": np.stack(
                        [
                            out_y1_2deg_agri[:, :, is_agri_is_past].squeeze(axis=2),
                            np.sum(out_y1_2deg_agri[:, :, ~is_agri_is_past], axis=2),
                            out_y1_2deg_ntrl,
                            out_y1_2deg_bare,
                        ],
                        axis=2,
                    )
                    / np.repeat(land_area_2deg_yx[:, :, None], 4, axis=2),
                }
                savemat(
                    os.path.join(
                        harm_dir_year, f"LandCoverFract.base{cfg.base_year}.2deg.mat"
                    ),
                    {"out_y1": out_lcf_2},
                )
                if not cfg.combine_crops:
                    agri_cf_2 = out_y1_2deg_agri[:, :, ~is_agri_is_past]
                    denom_2deg = np.sum(agri_cf_2, axis=2)
                    out_cf_2 = {
                        "varNames": np.array(list(lpjg_crops), dtype=object),
                        "maps_YXv": agri_cf_2
                        / np.repeat(denom_2deg[:, :, None], len(lpjg_crops), axis=2),
                    }
                    out_cf_2["maps_YXv"][np.isnan(out_cf_2["maps_YXv"])] = 0
                    out_cf_2["maps_YXv"][
                        np.repeat(denom_2deg[:, :, None] == 0, len(lpjg_crops), axis=2)
                    ] = 0
                    savemat(
                        os.path.join(harm_dir_year, f"CropFract.base{cfg.base_year}.2deg.mat"),
                        {"out_y1": out_cf_2},
                    )
                    out_fert_2 = {
                        "varNames": np.array(list(lpjg_crops), dtype=object),
                        "maps_YXv": 1e4 * out_y1_2deg_nfert,
                    }
                    savemat(
                        os.path.join(harm_dir_year, f"Fert.base{cfg.base_year}.2deg.mat"),
                        {"out_y1": out_fert_2},
                    )
                    ir_2deg_save = out_y1_2deg_irrig.copy()
                    ir_2deg_save[
                        np.repeat(denom_2deg[:, :, None] == 0, len(lpjg_crops), axis=2)
                    ] = 0
                    out_ir_2 = {
                        "varNames": np.array(list(lpjg_crops), dtype=object),
                        "maps_YXv": ir_2deg_save,
                    }
                    savemat(
                        os.path.join(harm_dir_year, f"Irrig.base{cfg.base_year}.2deg.mat"),
                        {"out_y1": out_ir_2},
                    )
            if cfg.save_2deg_txt:
                out_lcf_txt = {
                    "varNames": ["PASTURE", "CROPLAND", "NATURAL", "BARREN"],
                    "maps_YXv": np.stack(
                        [
                            out_y1_2deg_agri[:, :, is_agri_is_past].squeeze(axis=2),
                            np.sum(out_y1_2deg_agri[:, :, ~is_agri_is_past], axis=2),
                            out_y1_2deg_ntrl,
                            out_y1_2deg_bare,
                        ],
                        axis=2,
                    )
                    / np.repeat(land_area_2deg_yx[:, :, None], 4, axis=2),
                }
                out_array, out_header = maps2table(out_lcf_txt, list2map_2deg)
                save_table(
                    out_header,
                    out_array,
                    os.path.join(harm_dir_year, f"LandCoverFract.base{cfg.base_year}.2deg.txt"),
                    outPrec=cfg.out_prec,
                    outWidth=cfg.out_width,
                    delimiter=cfg.delimiter,
                    overwrite=cfg.overwrite,
                    fancy=cfg.fancy,
                    progress_step_pct=20,
                    verbose=False,
                )
                if not cfg.combine_crops:
                    agri_cf_2 = out_y1_2deg_agri[:, :, ~is_agri_is_past]
                    denom_2deg = np.sum(agri_cf_2, axis=2)
                    out_cf_txt = {
                        "varNames": list(lpjg_crops),
                        "maps_YXv": agri_cf_2
                        / np.repeat(denom_2deg[:, :, None], len(lpjg_crops), axis=2),
                    }
                    out_cf_txt["maps_YXv"][np.isnan(out_cf_txt["maps_YXv"])] = 0
                    out_array, out_header = maps2table(out_cf_txt, list2map_2deg)
                    save_table(
                        out_header,
                        out_array,
                        os.path.join(harm_dir_year, f"CropFract.base{cfg.base_year}.2deg.txt"),
                        outPrec=cfg.out_prec,
                        outWidth=cfg.out_width,
                        delimiter=cfg.delimiter,
                        overwrite=cfg.overwrite,
                        fancy=cfg.fancy,
                        progress_step_pct=20,
                        verbose=False,
                    )
                    out_fert_txt = {
                        "varNames": list(lpjg_crops),
                        "maps_YXv": 1e4 * out_y1_2deg_nfert,
                    }
                    out_array, out_header = maps2table(out_fert_txt, list2map_2deg)
                    save_table(
                        out_header,
                        out_array,
                        os.path.join(harm_dir_year, f"Fert.base{cfg.base_year}.2deg.txt"),
                        outPrec=cfg.out_prec,
                        outWidth=cfg.out_width,
                        delimiter=cfg.delimiter,
                        overwrite=cfg.overwrite,
                        fancy=cfg.fancy,
                        progress_step_pct=20,
                        verbose=False,
                    )
                    ir_2deg_txt_save = out_y1_2deg_irrig.copy()
                    ir_2deg_txt_save[
                        np.repeat(denom_2deg[:, :, None] == 0, len(lpjg_crops), axis=2)
                    ] = 0
                    out_ir_txt = {"varNames": list(lpjg_crops), "maps_YXv": ir_2deg_txt_save}
                    out_array, out_header = maps2table(out_ir_txt, list2map_2deg)
                    save_table(
                        out_header,
                        out_array,
                        os.path.join(harm_dir_year, f"Irrig.base{cfg.base_year}.2deg.txt"),
                        outPrec=cfg.out_prec,
                        outWidth=cfg.out_width,
                        delimiter=cfg.delimiter,
                        overwrite=cfg.overwrite,
                        fancy=cfg.fancy,
                        progress_step_pct=20,
                        verbose=False,
                    )

            print(f"  Done ({_toc_hms(time.perf_counter() - tic)}).")

            out_y0_agri_yxv = out_y1_agri
            out_y0_vegd_yx = np.sum(out_y1_agri, axis=2) + out_y1_ntrl
            out_y0_2deg_agri_yxv = out_y1_2deg_agri
            out_y0_2deg_vegd_yx = np.sum(out_y1_2deg_agri, axis=2) + out_y1_2deg_ntrl
            out_y0 = {
                "varNames": np.array(list(lu_names), dtype=object),
                "maps_YXv": np.concatenate(
                    [out_y1_agri, out_y1_ntrl[:, :, None], out_y1_bare[:, :, None]],
                    axis=2,
                ),
            }
            out_y0_2deg = {
                "varNames": np.array(list(lu_names), dtype=object),
                "maps_YXv": np.concatenate(
                    [
                        out_y1_2deg_agri,
                        out_y1_2deg_ntrl[:, :, None],
                        out_y1_2deg_bare[:, :, None],
                    ],
                    axis=2,
                ),
            }
            if not cfg.combine_crops:
                out_y0_2deg_nfert = out_y1_2deg_nfert
                out_y0_2deg_irrig = out_y1_2deg_irrig
                out_y0_nfert = out_y1_nfert
                out_y0_irrig = out_y1_irrig
                max_nfert_y0 = max_nfert_y1
                max_orig_nfert_y0 = max_orig_nfert_y1

            in_y0 = in_y1
            in_y0_2deg = in_y1_2deg
            in_y0_agri = in_y1_agri
            in_y0_2deg_agri = in_y1_2deg_agri
            in_y0_nfert = in_y1_nfert
            in_y0_irrig = in_y1_irrig
            in_y0_2deg_nfert = in_y1_2deg_nfert
            in_y0_2deg_irrig = in_y1_2deg_irrig
            bare_frac_y0 = (
                in_y1["maps_YXv"][:, :, lu_names.index("BARREN")] / land_area_yx
            )
            harm_dir_year = os.path.join(harm_dir, str(this_year))
            os.makedirs(harm_dir_year, exist_ok=True)
            post_path = os.path.join(harm_dir_year, f"post.base{cfg.base_year}.mat")
            in_y0_2deg_vegd_yx = np.sum(
                in_y0_2deg["maps_YXv"][:, :, not_bare], axis=2
            )
            out_y0_2deg_nfert_struct = (
                {"varNames": list(lpjg_crops), "maps_YXv": out_y0_2deg_nfert}
                if not cfg.combine_crops
                else None
            )
            out_y0_2deg_irrig_struct = (
                {"varNames": list(lpjg_crops), "maps_YXv": out_y0_2deg_irrig}
                if not cfg.combine_crops
                else None
            )
            post_vars = {
                "bareFrac_y0_YX": bare_frac_y0,
                "in_y0": in_y0,
                "in_y0_2deg": in_y0_2deg,
                "in_y0_2deg_agri_YXv": in_y0_2deg_agri,
                "in_y0_2deg_irrig": in_y0_2deg_irrig,
                "in_y0_2deg_nfert": in_y0_2deg_nfert,
                "in_y0_2deg_vegd_YX": in_y0_2deg_vegd_yx,
                "in_y0_agri_YXv": in_y0_agri,
                "in_y0_irrig": in_y0_irrig,
                "in_y0_nfert": in_y0_nfert,
                "latestPLUMin_2deg_irrig_YXv": latest_irrig,
                "latestPLUMin_2deg_nfert_YXv": latest_nfert,
                "max_harm_nfert_y0": max_harm_nfert_y0,
                "max_nfert_y0": max_nfert_y0 if not cfg.combine_crops else None,
                "max_orig_nfert_y0": max_orig_nfert_y0 if not cfg.combine_crops else None,
                "out_y0": out_y0,
                "out_y0_2deg": out_y0_2deg,
                "out_y0_2deg_agri_YXv": out_y0_2deg_agri_yxv,
                "out_y0_2deg_irrig": out_y0_2deg_irrig_struct,
                "out_y0_2deg_irrig_YXv": out_y0_2deg_irrig,
                "out_y0_2deg_nfert": out_y0_2deg_nfert_struct,
                "out_y0_2deg_nfert_YXv": out_y0_2deg_nfert,
                "out_y0_2deg_ntrl_YX": out_y0_2deg_ntrl_yx,
                "out_y0_2deg_vegd_YX": out_y0_2deg_vegd_yx,
                "out_y0_agri_YXv": out_y0_agri_yxv,
                "out_y0_irrig_YXv": out_y0_irrig,
                "out_y0_nfert_YXv": out_y0_nfert,
                "out_y0_vegd_YX": out_y0_vegd_yx,
            }
            _save_post_mat(post_path, post_vars)

            # Accumulate time-series for inline figures
            y_idx = years.index(this_year)
            if is_first_year_in_loop:
                year_list_ts = [years[0] - 1] + list(years)
                n_lu = len(lu_names)
                n_crops = len(lpjg_crops) if not cfg.combine_crops else 0
                ts_in = {
                    "yearList": year_list_ts,
                    "luNames": list(lu_names),
                    "cropNames": list(lpjg_crops) if not cfg.combine_crops else [],
                    "area_vy": np.full((n_lu, n_years + 1), np.nan),
                    "nfert_vy": np.full((n_crops, n_years + 1), np.nan),
                    "irrig_vy": np.full((n_crops, n_years + 1), np.nan),
                }
                ts_out = {k: (v.copy() if isinstance(v, np.ndarray) else list(v) if isinstance(v, list) else v) for k, v in ts_in.items()}
                _save_ts(ts_in, in_y0, in_y0_nfert["maps_YXv"] if not cfg.combine_crops else None, in_y0_irrig["maps_YXv"] if not cfg.combine_crops else None, 0, lu_names, lpjg_crops)
                _save_ts(ts_out, out_y0, out_y0_nfert if not cfg.combine_crops else None, out_y0_irrig if not cfg.combine_crops else None, 0, lu_names, lpjg_crops)

            _save_ts(ts_in, in_y1, in_y1_nfert["maps_YXv"] if not cfg.combine_crops else None, in_y1_irrig["maps_YXv"] if not cfg.combine_crops else None, y_idx + 1, lu_names, lpjg_crops)
            out_y1_full = {
                "varNames": list(lu_names),
                "maps_YXv": np.concatenate([out_y1_agri, out_y1_ntrl[:, :, None], out_y1_bare[:, :, None]], axis=2),
            }
            _save_ts(ts_out, out_y1_full, out_y1_nfert if not cfg.combine_crops else None, out_y1_irrig if not cfg.combine_crops else None, y_idx + 1, lu_names, lpjg_crops)

            is_first_year_in_loop = False

        # Generate inline time-series figures
        if ts_in.get("area_vy") is not None:
            fig_dir = harm_dir + "_figs"
            os.makedirs(fig_dir, exist_ok=True)
            _make_inline_figs(ts_in, ts_out, lu_names, lpjg_crops if not cfg.combine_crops else [], is_crop, fig_dir)

    print("Done")


def _get_set_area_to_zero(
    lu_names: Sequence[str],
    in_y0_yxv: np.ndarray,
    in_y1_yxv: np.ndarray,
    out_y0_yxv: np.ndarray,
    this_year: int,
    base_year: int,
    conserv_tol_area: float,
) -> tuple[np.ndarray, bool]:
    set_area_to_zero = np.zeros(len(lu_names), dtype=bool)
    updated_in_y0 = False
    for i, this_lu in enumerate(lu_names):
        out0 = out_y0_yxv[:, :, i]
        in0 = in_y0_yxv[:, :, i]
        in1 = in_y1_yxv[:, :, i]
        in_diff = in1 - in0
        glob_loss = np.sum(in_diff[in_diff < 0])
        glob_area = np.sum(out0)
        net_glob_chg = np.sum(in_diff)
        if -glob_loss > glob_area:
            if this_lu == "Miscanthus" and this_year - 1 == base_year:
                if np.sum(in1[in0 > 0]) == 0:
                    print(
                        "Warning: PLUM has Miscanthus in 2010, and some gridcells lose Miscanthus from 2010 "
                        "to 2011. This cannot be satisfied using my 2010 maps because I have no Miscanthus. "
                        "However, it shouldn't be a problem because the grid cells PLUM specifies as losing "
                        "Miscanthus lose ALL it."
                    )
                else:
                    diff = np.sum(in1[in0 > 0])
                    print(
                        "Warning: PLUM has Miscanthus in 2010, and some gridcells lose Miscanthus from 2010 "
                        "to 2011. This cannot be satisfied using my 2010 maps because I have no Miscanthus. "
                        "But not all grid cells PLUM specifies as losing Miscanthus lose ALL of it! "
                        f"Difference of {diff:0.1g} m2. Ignoring."
                    )
                tmp = in_y0_yxv[:, :, i].copy()
                tmp[in_diff < 0] = in1[in_diff < 0]
                in_y0_yxv[:, :, i] = tmp
                updated_in_y0 = True
            elif 0 <= net_glob_chg < conserv_tol_area:
                print(
                    f"Warning: Total global loss of {this_lu} ({-glob_loss:0.4e} m2) exceeds its y0 area "
                    f"({glob_area:0.4e} m2). LIKELY WILL cause infinite loop in ringRedist.\n"
                )
            elif glob_area + np.sum(in_diff[in_diff > 0]) > -glob_loss:
                glob_area_plus_gain = glob_area + np.sum(in_diff[in_diff > 0])
                print(
                    f"Warning: Total global loss of {this_lu} ({-glob_loss:0.4e} m2) exceeds its y0 area "
                    f"({glob_area:0.4e} m2), but not y0 area + total global gain ({glob_area_plus_gain:0.4e} m2). "
                    "(((May))) cause infinite loop in ringRedist.\n"
                )
            else:
                glob_area_plus_gain = glob_area + np.sum(in_diff[in_diff > 0])
                print(
                    f"Warning: Total global loss of {this_lu} ({-glob_loss:0.4e} m2) exceeds its y0 area "
                    f"({glob_area:0.4e} m2) AND ALSO y0 area + total global gain ({glob_area_plus_gain:0.4e} m2). "
                    "Setting area to zero.\n"
                )
                set_area_to_zero[i] = True
    return set_area_to_zero, updated_in_y0


def _check_discrepancy(tmp_in: float, tmp_out: float, this_var: str, this_name: str) -> None:
    err_pct = (tmp_in - tmp_out) / tmp_out * 100 if tmp_out != 0 else np.nan
    err_sign = "+" if err_pct > 0 else ""
    print(f"        {this_var} {this_name}: {err_sign}{err_pct:0.1f}%")


def _save_post_mat(path: str, post_vars: dict) -> None:
    struct_vars = {
        "in_y0",
        "in_y0_2deg",
        "in_y0_nfert",
        "in_y0_irrig",
        "in_y0_2deg_nfert",
        "in_y0_2deg_irrig",
        "out_y0",
        "out_y0_2deg",
        "out_y0_2deg_nfert",
        "out_y0_2deg_irrig",
    }
    with h5py.File(path, "w") as f:
        for key, val in post_vars.items():
            if val is None:
                continue
            if key in struct_vars:
                _write_struct(f, key, val)
            else:
                arr = _coerce_array(val)
                if key.startswith("max_"):
                    arr = np.asarray(arr).reshape(1, -1)
                    _write_dataset(f, key, arr, matlab_order=False)
                else:
                    _write_dataset(f, key, arr, matlab_order=True)


def _find_post_restart_file(harm_dir: str, prev_year: int, base_year: int) -> str | None:
    candidate_year_dir = os.path.join(harm_dir, str(prev_year))
    cand_1 = os.path.join(candidate_year_dir, f"post.base{base_year}.mat")
    cand_2 = os.path.join(harm_dir, f"{prev_year}post.base{base_year}.mat")
    if os.path.isfile(cand_1):
        return cand_1
    if os.path.isfile(cand_2):
        return cand_2
    return None


def _load_post_mat(path: str) -> dict:
    post: dict = {}
    with h5py.File(path, "r") as f:
        for key in f.keys():
            if key == "#refs#":
                continue
            obj = f[key]
            if isinstance(obj, h5py.Group):
                post[key] = _read_struct(obj)
            else:
                post[key] = _read_dataset(obj, key)
    post = _rename_post_keys(post)
    return post


def _rename_post_keys(post: dict) -> dict:
    renamed = dict(post)
    for key in list(post.keys()):
        if "nfert_2deg" in key:
            new_key = key.replace("nfert_2deg", "2deg_nfert")
        elif "irrig_2deg" in key:
            new_key = key.replace("irrig_2deg", "2deg_irrig")
        else:
            continue
        if new_key not in renamed:
            renamed[new_key] = renamed[key]
        renamed.pop(key, None)
    return renamed


def _read_struct(group: h5py.Group) -> dict:
    maps = group.get("maps_YXv", None)
    var_names = group.get("varNames", None)
    out: dict = {}
    if maps is not None:
        out["maps_YXv"] = _from_matlab_order(maps[...])
    if var_names is not None:
        out["varNames"] = _read_varnames(var_names)
    return out


def _read_varnames(ds: h5py.Dataset) -> list[str]:
    ref_dtype = h5py.check_dtype(ref=ds.dtype)
    if ref_dtype is not None:
        names = []
        for i in range(ds.shape[0]):
            ref = ds[i, 0]
            if not isinstance(ref, h5py.Reference):
                names.append(str(ref))
                continue
            codes = ds.file[ref][...]
            names.append(_decode_char_codes(codes))
        return names
    data = ds[...]
    if data.dtype.kind in {"S", "U"}:
        return [str(x).strip() for x in data.ravel()]
    return [str(x) for x in data.ravel()]


def _decode_char_codes(codes: np.ndarray) -> str:
    flat = np.asarray(codes).ravel()
    chars = [chr(int(c)) for c in flat if int(c) != 0]
    return "".join(chars)


def _read_dataset(ds: h5py.Dataset, key: str) -> np.ndarray:
    data = ds[...]
    if key.startswith("max_"):
        return np.asarray(data).squeeze()
    return _from_matlab_order(np.asarray(data))


def _from_matlab_order(arr: np.ndarray) -> np.ndarray:
    if arr.ndim == 2:
        return arr.T
    if arr.ndim == 3:
        return arr.transpose(2, 1, 0)
    return arr


def _get_post_maps(post: dict, struct_key: str, yxv_key: str) -> np.ndarray | None:
    if yxv_key in post and isinstance(post[yxv_key], np.ndarray):
        return post[yxv_key]
    struct = post.get(struct_key, None)
    if isinstance(struct, dict):
        return struct.get("maps_YXv", None)
    return None


def _write_struct(f: h5py.File, name: str, data: dict) -> None:
    if data is None:
        return
    var_names = data.get("varNames", [])
    maps = data.get("maps_YXv", None)
    if maps is None:
        return
    grp = f.create_group(name)
    _write_dataset(grp, "maps_YXv", maps, matlab_order=True)
    _write_varnames(grp, var_names)


def _write_varnames(group: h5py.Group, var_names: Sequence[str]) -> None:
    ref_dtype = h5py.special_dtype(ref=h5py.Reference)
    ds = group.create_dataset("varNames", shape=(len(var_names), 1), dtype=ref_dtype)
    refs = group.file.require_group("#refs#")
    for i, name in enumerate(var_names):
        codes = np.array([[ord(c)] for c in str(name)], dtype=np.uint16)
        dset = refs.require_dataset(
            f"varName_{i}", shape=codes.shape, dtype=codes.dtype
        )
        dset[...] = codes
        ds[i, 0] = dset.ref


def _write_dataset(
    parent: h5py.File | h5py.Group,
    name: str,
    arr: np.ndarray,
    matlab_order: bool,
) -> None:
    data = _coerce_array(arr)
    if matlab_order:
        data = _to_matlab_order(data)
    parent.create_dataset(name, data=data)


def _to_matlab_order(arr: np.ndarray) -> np.ndarray:
    if arr.ndim == 2:
        return arr.T
    if arr.ndim == 3:
        return arr.transpose(2, 1, 0)
    return arr


def _coerce_array(val: object) -> np.ndarray:
    if isinstance(val, np.ndarray):
        return val
    return np.asarray(val)


def _debug_global_areas(
    maps_a_yxv: np.ndarray,
    maps_b_yxv: np.ndarray,
    msg_intro: str,
    msg_a: str,
    msg_b: str,
    lu_names: Sequence[str],
    is_crop: np.ndarray,
    is_agri: np.ndarray,
    db_crop: int | str | None,
    this_year: int,
) -> None:
    msg_width = 45
    units = "m2"
    print("")
    print(msg_intro)
    if db_crop not in (None, ""):
        if isinstance(db_crop, str):
            db_crop = lu_names.index(db_crop)
        tmp0 = np.sum(maps_a_yxv[:, :, db_crop])
        tmp1 = np.sum(maps_b_yxv[:, :, db_crop])
        print(f"{msg_a}_{this_year-1} glob {lu_names[db_crop]} area ({units}):".ljust(msg_width) + f" {tmp0:0.4e}")
        print(f"{msg_b}_{this_year} glob {lu_names[db_crop]} area ({units}):".ljust(msg_width) + f" {tmp1:0.4e}")
        print(f"diff ({units}):".ljust(msg_width) + f" {tmp1-tmp0:0.4e}\n")
    tmp0 = np.sum(maps_a_yxv[:, :, is_crop])
    tmp1 = np.sum(maps_b_yxv[:, :, is_crop])
    print(f"{msg_a}_{this_year-1} glob crop area ({units}):".ljust(msg_width) + f" {tmp0:0.4e}")
    print(f"{msg_b}_{this_year} glob crop area ({units}):".ljust(msg_width) + f" {tmp1:0.4e}")
    print(f"diff ({units}):".ljust(msg_width) + f" {tmp1-tmp0:0.4e}\n")


def _debug_global_mgmts(
    maps_a_mgmt_yxv: np.ndarray,
    maps_b_mgmt_yxv: np.ndarray,
    unmet_mgmt_yxv: np.ndarray,
    maps_a_area_yxv: np.ndarray,
    maps_b_area_yxv: np.ndarray,
    units: str,
    convfact: float,
    msg_intro: str,
    msg_a: str,
    msg_b: str,
    msg_mgmt: str,
    lpjg_crops: Sequence[str],
    db_crop: int | str | None,
    this_year: int,
) -> None:
    msg_width = 45
    maps_a_yxv = convfact * (maps_a_mgmt_yxv * maps_a_area_yxv)
    maps_b_yxv = convfact * (maps_b_mgmt_yxv * maps_b_area_yxv + unmet_mgmt_yxv)
    print("")
    print(msg_intro)
    if db_crop not in (None, ""):
        if isinstance(db_crop, str):
            db_crop = lpjg_crops.index(db_crop)
        tmp0 = np.sum(maps_a_yxv[:, :, db_crop])
        tmp1 = np.sum(maps_b_yxv[:, :, db_crop])
        print(f"{msg_a}_{this_year-1} glob {lpjg_crops[db_crop]} {msg_mgmt} ({units}):".ljust(msg_width) + f" {tmp0:0.4e}")
        print(f"{msg_b}_{this_year} glob {lpjg_crops[db_crop]} {msg_mgmt} ({units}):".ljust(msg_width) + f" {tmp1:0.4e}")
        print(f"diff ({units}):".ljust(msg_width) + f" {tmp1-tmp0:0.4e}\n")
    tmp0 = np.sum(maps_a_yxv)
    tmp1 = np.sum(maps_b_yxv)
    print(f"{msg_a}_{this_year-1} glob total {msg_mgmt} ({units}):".ljust(msg_width) + f" {tmp0:0.4e}")
    print(f"{msg_b}_{this_year} glob total {msg_mgmt} ({units}):".ljust(msg_width) + f" {tmp1:0.4e}")
    print(f"diff ({units}):".ljust(msg_width) + f" {tmp1-tmp0:0.4e}\n")


def _check_preserved_global_deltas(
    out_y0_a_yxv: np.ndarray,
    out_y1_a_yxv: np.ndarray,
    out_y0_b_yxv: np.ndarray,
    out_y1_b_yxv: np.ndarray,
    conserv_tol_pct: float,
    var_names: Sequence[str],
    check_name: str,
    area_type: str,
) -> None:
    delta_2deg = np.sum(out_y1_a_yxv - out_y0_a_yxv)
    delta = np.sum(out_y1_b_yxv - out_y0_b_yxv)
    if delta_2deg != 0:
        delta_diff = delta - delta_2deg
        delta_diff_pct = 100 * delta_diff / abs(delta_2deg)
        if abs(delta_diff_pct) > conserv_tol_pct:
            print(
                f"Warning: {check_name} changes ∆ {area_type} area: "
                f"Diff {delta_diff} ({delta_diff_pct:0.2f}%, tolerance {conserv_tol_pct:0.2f}%)\n"
            )
    for v, name in enumerate(var_names):
        delta_2deg = np.sum(out_y1_a_yxv[:, :, v] - out_y0_a_yxv[:, :, v])
        delta = np.sum(out_y1_b_yxv[:, :, v] - out_y0_b_yxv[:, :, v])
        if delta_2deg == 0:
            continue
        delta_diff = delta - delta_2deg
        delta_diff_pct = 100 * delta_diff / abs(delta_2deg)
        if abs(delta_diff_pct) > conserv_tol_pct:
            print(
                f"Warning: {check_name} changes ∆ {name} area: "
                f"Diff {delta_diff} ({delta_diff_pct:0.2f}%, tolerance {conserv_tol_pct:0.2f}%)\n"
            )


def _align_maps_by_names(
    maps_yxv: np.ndarray, var_names: Sequence[str], target_names: Sequence[str]
) -> np.ndarray:
    out = np.zeros((maps_yxv.shape[0], maps_yxv.shape[1], len(target_names)), dtype=maps_yxv.dtype)
    name_to_idx = {v: i for i, v in enumerate(var_names)}
    for i, name in enumerate(target_names):
        if name in name_to_idx:
            out[:, :, i] = maps_yxv[:, :, name_to_idx[name]]
    return out


def _read_min_natural_rate(plum_config_file: str) -> float:
    min_natural_rate = None
    with open(plum_config_file, "r") as handle:
        for line in handle:
            if line.strip().startswith("MIN_NATURAL_RATE"):
                parts = line.strip().split("=", 1)
                if len(parts) == 2:
                    min_natural_rate = parts[1].strip()
                break
    if min_natural_rate is None or min_natural_rate == "":
        print(f"Warning: MIN_NATURAL_RATE not found in {plum_config_file}; using 0")
        return 0.0
    try:
        return float(min_natural_rate)
    except ValueError:
        print(f"Warning: MIN_NATURAL_RATE in {plum_config_file} is not numeric; using 0")
        return 0.0


def _toc_hms(toc_sec: float) -> str:
    if toc_sec >= 3600:
        hours = int(toc_sec // 3600)
        minutes = int((toc_sec % 3600) // 60)
        seconds = int(toc_sec % 60)
        return f"{hours:02d}:{minutes:02d}:{seconds:02d}"
    if toc_sec >= 60:
        minutes = int(toc_sec // 60)
        seconds = int(toc_sec % 60)
        return f"{minutes:02d}:{seconds:02d}"
    return f"{round(toc_sec * 1e1) * 1e-1} seconds"


def _normalize_dirs(dirs: Sequence[str], base: str | None) -> List[str]:
    out = []
    for d in dirs:
        if base and not os.path.isabs(d):
            out.append(os.path.join(base, d))
        else:
            out.append(d)
    return out


def _remove_trailing_slash(path: str) -> str:
    while path.endswith(os.sep):
        path = path[:-1]
    return path


def _get_harm_dir(in_dir: str, fruitveg_sugar_2oil: bool, combine_crops: bool) -> str:
    harm_dir = in_dir
    if fruitveg_sugar_2oil:
        harm_dir = f"{harm_dir}.fvs2oil"
    if combine_crops:
        harm_dir = f"{harm_dir}.combineCrops"
    return harm_dir


def _get_harm_dirs(plum_dirs: Sequence[str], fruitveg_sugar_2oil: bool, combine_crops: bool) -> List[str]:
    return [
        _get_harm_dir(f"{d}.harm", fruitveg_sugar_2oil, combine_crops) for d in plum_dirs
    ]


def _check_dirs(dir_list: Sequence[str], rw: str) -> List[str]:
    checked: List[str] = []
    for d in dir_list:
        this_dir = _remove_trailing_slash(d)
        if not os.path.isdir(this_dir):
            if "w" not in rw:
                raise RuntimeError(f"{this_dir} not found")
            os.makedirs(this_dir, exist_ok=True)
        this_dir = os.path.abspath(this_dir)
        if "r" in rw and not os.access(this_dir, os.R_OK):
            raise RuntimeError(f"{this_dir} is not readable!")
        if "w" in rw and not os.access(this_dir, os.W_OK):
            raise RuntimeError(f"{this_dir} is not writeable!")
        if not os.path.isdir(this_dir):
            raise RuntimeError(f"{this_dir} is not a directory!")
        checked.append(this_dir)
    return checked


def _save_ts(
    ts: dict,
    s_area: dict,
    nfert_yxv: np.ndarray | None,
    irrig_yxv: np.ndarray | None,
    y: int,
    lu_names: Sequence[str],
    crop_names: Sequence[str],
) -> None:
    """Accumulate global totals into time-series struct (matches MATLAB save_to_timeseries_struct)."""
    ts["area_vy"][:, y] = np.nansum(np.nansum(s_area["maps_YXv"], axis=0), axis=0)
    if nfert_yxv is not None and irrig_yxv is not None:
        var_names = list(s_area.get("varNames", lu_names))
        crop_idx_in_lu = [var_names.index(c) for c in crop_names if c in var_names]
        crop_areas = s_area["maps_YXv"][:, :, crop_idx_in_lu]
        ts["nfert_vy"][:, y] = np.nansum(np.nansum(nfert_yxv * crop_areas, axis=0), axis=0)
        ts["irrig_vy"][:, y] = np.nansum(np.nansum(irrig_yxv * crop_areas, axis=0), axis=0)


def _make_inline_figs(
    ts_in: dict,
    ts_out: dict,
    lu_names: Sequence[str],
    crop_names: Sequence[str],
    is_crop: np.ndarray,
    fig_dir: str,
) -> None:
    """Generate inline time-series figures (matches MATLAB PLUMharm.m lines 1376-1423)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available; skipping inline figures.")
        return

    year_list = ts_in["yearList"]
    legend_ts = ["From PLUM", "Harmonized"]

    def _make_fig(ts_orig, ts_harm, names, units, title_word, file_word):
        n_vars = len(names)
        ny = int(np.ceil(n_vars / 2))
        fig, axes = plt.subplots(ny, 2, figsize=(15, 4 * ny), constrained_layout=True)
        axes = np.atleast_2d(axes)
        for v in range(n_vars):
            ax = axes[v // 2, v % 2]
            ax.plot(year_list, ts_orig[v, :], "--", linewidth=1, label=legend_ts[0])
            ax.plot(year_list, ts_harm[v, :], "-", linewidth=1, label=legend_ts[1])
            ax.set_title(f"{title_word}: {names[v]}")
            ax.set_ylabel(units)
            ax.legend(loc="best", fontsize=8)
        for v in range(n_vars, ny * 2):
            axes[v // 2, v % 2].set_visible(False)
        fig.savefig(os.path.join(fig_dir, f"timeSeries_{file_word}.pdf"))
        plt.close(fig)

    # LU areas
    lu_name_list = list(lu_names)
    not_crop = ~is_crop
    tmp_names = ["CROPLAND"] + [lu_name_list[i] for i in range(len(lu_name_list)) if not_crop[i]]
    n_tmp = len(tmp_names)
    ts_orig_lu = np.full((n_tmp, len(year_list)), np.nan)
    ts_harm_lu = np.full((n_tmp, len(year_list)), np.nan)
    ts_orig_lu[0, :] = np.sum(ts_in["area_vy"][is_crop, :], axis=0)
    ts_orig_lu[1:, :] = ts_in["area_vy"][not_crop, :]
    ts_harm_lu[0, :] = np.sum(ts_out["area_vy"][is_crop, :], axis=0)
    ts_harm_lu[1:, :] = ts_out["area_vy"][not_crop, :]
    _make_fig(ts_orig_lu * 1e-12, ts_harm_lu * 1e-12, tmp_names, "Million km²", "Area", "landUse")

    if not crop_names:
        return

    # Crop areas
    crop_area_in = ts_in["area_vy"][is_crop, :] * 1e-12
    crop_area_out = ts_out["area_vy"][is_crop, :] * 1e-12
    _make_fig(crop_area_in, crop_area_out, list(crop_names), "Million km²", "Area", "crops")

    # Fertilization
    nfert_in = ts_in["nfert_vy"] * 1e-3 * 1e-6
    nfert_out = ts_out["nfert_vy"] * 1e-3 * 1e-6
    _make_fig(nfert_in, nfert_out, list(crop_names), "Mt N", "Fert.", "nfert")

    # Irrigation
    _make_fig(ts_in["irrig_vy"], ts_out["irrig_vy"], list(crop_names), "intensity × area", "Irrigation", "irrig")
