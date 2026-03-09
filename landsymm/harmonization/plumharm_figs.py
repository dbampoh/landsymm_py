"""Diagnostic figure generation (partial)."""
from __future__ import annotations

import os
import shutil
from typing import Any, Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .plumharm_figs_options import PlumHarmFigsConfig
from .plumharm_import_ref import import_ref_data
from .plumharm_pp_read_plum import pp_read_plum


def run_plumharm_figs(cfg: PlumHarmFigsConfig) -> Dict[str, np.ndarray]:
    """Prepare data arrays for PLUMharm diagnostic figures."""
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
    harm_dirs = _check_dirs(harm_dirs, "r")

    if not os.path.isdir(cfg.harms_figs_dir):
        os.makedirs(cfg.harms_figs_dir, exist_ok=True)
    harms_figs_dir = os.path.abspath(cfg.harms_figs_dir)

    year_list_harm = list(range(cfg.year1, cfg.yearN + 1))
    year_list_orig = [year_list_harm[0] - 1] + year_list_harm

    if cfg.year_list_baseline_to_plot is None:
        year_list_baseline = year_list_harm
    else:
        year_list_baseline = list(cfg.year_list_baseline_to_plot)

    if cfg.this_ver not in {"", "2deg.", "orig."}:
        raise RuntimeError(f"Unrecognized value of this_ver: {cfg.this_ver}")

    if cfg.three_years is None:
        three_years = [year_list_harm[0], year_list_harm[len(year_list_harm) // 2], year_list_harm[-1]]
    else:
        three_years = list(cfg.three_years)
    if min(three_years) < min(year_list_harm) or max(three_years) > max(year_list_harm):
        raise RuntimeError("three_years must be entirely within year_list_harm")

    ref = import_ref_data(
        _to_plumharm_config(cfg),
        do_harm=False,
        year_list_baseline_lu_to_plot=year_list_baseline,
    )
    ref["geodata_dir"] = cfg.geodata_dir
    ref["template_dir"] = os.path.join(cfg.this_dir or os.getcwd(), "data", "templates")

    mask_yx = ref["mask_YX"]
    n_lu = len(ref["lu_names"])
    n_runs = len(plum_dirs)

    is2deg = cfg.this_ver == "2deg."
    if is2deg:
        ny, nx = 90, 180
    else:
        ny, nx = 360, 720

    n_cells = np.sum(~mask_yx)
    plum_orig = np.full((n_cells, n_lu, len(year_list_orig), n_runs), np.nan, dtype=np.float32)
    plum_harm = np.full((n_cells, n_lu, len(year_list_harm), n_runs), np.nan, dtype=np.float32)
    if not cfg.combine_crops:
        n_crops = len(ref["lpjg_crops"])
        plum_orig_nfert = np.full((n_cells, n_crops, len(year_list_orig), n_runs), np.nan, dtype=np.float32)
        plum_orig_irrig = np.full((n_cells, n_crops, len(year_list_orig), n_runs), np.nan, dtype=np.float32)
        plum_harm_nfert = np.full((n_cells, n_crops, len(year_list_harm), n_runs), np.nan, dtype=np.float32)
        plum_harm_irrig = np.full((n_cells, n_crops, len(year_list_harm), n_runs), np.nan, dtype=np.float32)
    else:
        plum_orig_nfert = plum_orig_irrig = plum_harm_nfert = plum_harm_irrig = None

    for r, (plum_dir, harm_dir) in enumerate(zip(plum_dirs, harm_dirs)):
        if cfg.combine_crops:
            s_out, _, _ = pp_read_plum(
                plum_dir,
                cfg.base_year,
                year_list_orig,
                ref["landArea_YX"] if not is2deg else ref["landArea_2deg_YX"],
                ref["lu_names"],
                ref["plum_to_lpjg"],
                ref["lpjg_crops"],
                is2deg,
                [],
                cfg.norm2extra,
                None,
                cfg.this_ver,
                True,
                cfg.fruitveg_sugar_2oil,
                cfg.allow_unveg,
            )
        else:
            s_out, s_nf, s_ir = pp_read_plum(
                plum_dir,
                cfg.base_year,
                year_list_orig,
                ref["landArea_YX"] if not is2deg else ref["landArea_2deg_YX"],
                ref["lu_names"],
                ref["plum_to_lpjg"],
                ref["lpjg_crops"],
                is2deg,
                [],
                cfg.norm2extra,
                None,
                cfg.this_ver,
                True,
                cfg.fruitveg_sugar_2oil,
                cfg.allow_unveg,
            )
        year_indices = [list(s_out["yearList"]).index(y) for y in year_list_orig]
        s_out["maps_YXvy"] = s_out["maps_YXvy"][:, :, :, year_indices]
        incl = np.repeat(~mask_yx[:, :, None], s_out["maps_YXvy"].shape[2], axis=2)
        plum_orig[:, :, :, r] = s_out["maps_YXvy"][incl].reshape(n_cells, n_lu, len(year_list_orig))

        if not cfg.combine_crops:
            s_nf["maps_YXvy"] = s_nf["maps_YXvy"][:, :, :, year_indices]
            s_ir["maps_YXvy"] = s_ir["maps_YXvy"][:, :, :, year_indices]
            incl_nf = np.repeat(~mask_yx[:, :, None], s_nf["maps_YXvy"].shape[2], axis=2)
            plum_orig_nfert[:, :, :, r] = s_nf["maps_YXvy"][incl_nf].reshape(
                n_cells, n_crops, len(year_list_orig)
            )
            plum_orig_irrig[:, :, :, r] = s_ir["maps_YXvy"][incl_nf].reshape(
                n_cells, n_crops, len(year_list_orig)
            )

        if cfg.combine_crops:
            s_out, _, _ = pp_read_plum(
                harm_dir,
                cfg.base_year,
                year_list_harm,
                ref["landArea_YX"] if not is2deg else ref["landArea_2deg_YX"],
                ref["lu_names"],
                ref["plum_to_lpjg"],
                ref["lpjg_crops"],
                is2deg,
                [],
                0,
                None,
                cfg.this_ver,
                False,
                cfg.fruitveg_sugar_2oil,
                cfg.allow_unveg,
            )
        else:
            s_out, s_nf, s_ir = pp_read_plum(
                harm_dir,
                cfg.base_year,
                year_list_harm,
                ref["landArea_YX"] if not is2deg else ref["landArea_2deg_YX"],
                ref["lu_names"],
                ref["plum_to_lpjg"],
                ref["lpjg_crops"],
                is2deg,
                [],
                0,
                None,
                cfg.this_ver,
                False,
                cfg.fruitveg_sugar_2oil,
                cfg.allow_unveg,
            )
        year_indices = [list(s_out["yearList"]).index(y) for y in year_list_harm]
        s_out["maps_YXvy"] = s_out["maps_YXvy"][:, :, :, year_indices]
        incl = np.repeat(~mask_yx[:, :, None], s_out["maps_YXvy"].shape[2], axis=2)
        plum_harm[:, :, :, r] = s_out["maps_YXvy"][incl].reshape(n_cells, n_lu, len(year_list_harm))

        if not cfg.combine_crops:
            s_nf["maps_YXvy"] = s_nf["maps_YXvy"][:, :, :, year_indices]
            s_ir["maps_YXvy"] = s_ir["maps_YXvy"][:, :, :, year_indices]
            incl_nf = np.repeat(~mask_yx[:, :, None], s_nf["maps_YXvy"].shape[2], axis=2)
            plum_harm_nfert[:, :, :, r] = s_nf["maps_YXvy"][incl_nf].reshape(
                n_cells, n_crops, len(year_list_harm)
            )
            plum_harm_irrig[:, :, :, r] = s_ir["maps_YXvy"][incl_nf].reshape(
                n_cells, n_crops, len(year_list_harm)
            )

    if year_list_harm[0] == cfg.year1 and year_list_orig[0] == cfg.year1 - 1:
        plum_harm = np.concatenate([plum_orig[:, :, 0:1, :], plum_harm], axis=2)
        if not cfg.combine_crops:
            plum_harm_nfert = np.concatenate(
                [plum_orig_nfert[:, :, 0:1, :], plum_harm_nfert], axis=2
            )
            plum_harm_irrig = np.concatenate(
                [plum_orig_irrig[:, :, 0:1, :], plum_harm_irrig], axis=2
            )
        year_list_harm = [year_list_orig[0]] + year_list_harm

    if cfg.save_geotiffs:
        _write_geotiffs(
            harms_figs_dir,
            ref,
            plum_orig,
            plum_harm,
            year_list_orig,
            year_list_harm,
            cfg.this_ver,
            cfg.runlist_legend,
        )

    ts_orig_lu = np.nansum(plum_orig, axis=0)
    ts_harm_lu = np.nansum(plum_harm, axis=0)
    _make_crops_timeseries_fig(
        harms_figs_dir,
        ref["lu_names"],
        year_list_orig,
        year_list_harm,
        ts_orig_lu,
        ts_harm_lu,
        cfg.runlist_legend,
        suffix="landUse",
        title_word="Area",
        units="m\u00b2",
    )
    if not cfg.combine_crops:
        ts_orig_crops = np.nansum(plum_orig, axis=0)
        ts_harm_crops = np.nansum(plum_harm, axis=0)
        _make_crops_timeseries_fig(
            harms_figs_dir,
            ref["lpjg_crops"],
            year_list_orig,
            year_list_harm,
            ts_orig_crops,
            ts_harm_crops,
            cfg.runlist_legend,
            suffix="crops",
            title_word="Area",
            units="m\u00b2",
        )
        is_crop = ref["isCrop"]
        cf_kg2mt = 1e-3 * 1e-6
        ts_orig_nfert = cf_kg2mt * np.nansum(
            plum_orig[:, is_crop, :, :] * plum_orig_nfert, axis=0
        )
        ts_harm_nfert = cf_kg2mt * np.nansum(
            plum_harm[:, is_crop, :, :] * plum_harm_nfert, axis=0
        )
        _make_crops_timeseries_fig(
            harms_figs_dir,
            ref["lpjg_crops"],
            year_list_orig,
            year_list_harm,
            ts_orig_nfert,
            ts_harm_nfert,
            cfg.runlist_legend,
            suffix="nfert",
            title_word="Fert.",
            units="Mt N",
        )
        ts_orig_irrig_ts = np.nansum(
            plum_orig[:, is_crop, :, :] * plum_orig_irrig, axis=0
        )
        ts_harm_irrig_ts = np.nansum(
            plum_harm[:, is_crop, :, :] * plum_harm_irrig, axis=0
        )
        _make_crops_timeseries_fig(
            harms_figs_dir,
            ref["lpjg_crops"],
            year_list_orig,
            year_list_harm,
            ts_orig_irrig_ts,
            ts_harm_irrig_ts,
            cfg.runlist_legend,
            suffix="irrig",
            title_word="Irrigation",
            units="intensity \u00d7 area",
        )

    _plot_delta_maps(
        harms_figs_dir,
        ref["lu_names"],
        ref["mask_YX"],
        ref["list2map_2deg"] if is2deg else ref["list2map"],
        ny,
        nx,
        plum_orig,
        plum_harm,
        year_list_orig,
        year_list_harm,
    )
    _plot_three_year_maps(
        harms_figs_dir,
        ref,
        plum_orig,
        plum_harm,
        year_list_orig,
        year_list_harm,
        three_years,
        cfg.runlist_legend,
        is2deg,
        ny,
        nx,
    )
    _plot_change_maps(
        harms_figs_dir,
        ref,
        plum_orig,
        plum_harm,
        year_list_orig,
        year_list_harm,
        three_years,
        cfg.runlist_legend,
        is2deg,
        ny,
        nx,
    )
    _plot_diff_maps_three_years(
        harms_figs_dir,
        ref,
        plum_orig,
        plum_harm,
        year_list_orig,
        year_list_harm,
        three_years,
        cfg.runlist_legend,
        is2deg,
        ny,
        nx,
    )
    _plot_diff_maps_one_year_all_runs(
        harms_figs_dir,
        ref,
        plum_orig,
        plum_harm,
        year_list_orig,
        year_list_harm,
        three_years[1] if len(three_years) > 1 else three_years[0],
        cfg.runlist_legend,
        is2deg,
        ny,
        nx,
    )

    _plot_scatter_baseline_vs_orig(
        harms_figs_dir,
        ref,
        plum_orig,
        plum_harm,
        year_list_orig,
        year_list_harm,
        is2deg,
        ny,
        nx,
        runlist_legend=cfg.runlist_legend,
    )
    _plot_scatter_orig_vs_harm_delta(
        harms_figs_dir,
        ref,
        plum_orig,
        plum_harm,
        year_list_orig,
        year_list_harm,
        runlist_legend=cfg.runlist_legend,
    )

    _harm_by_numbers(
        harms_figs_dir,
        ref,
        plum_orig,
        plum_harm,
        year_list_orig,
        year_list_harm,
        cfg.runlist_legend,
        cfg.combine_crops,
    )

    _plot_harm_effect_nonagri(
        harms_figs_dir,
        ref,
        plum_orig,
        plum_harm,
        year_list_orig,
        cfg.runlist_legend,
        cfg.timeseries_legend_loc,
    )
    _plot_bundle_timeseries(
        harms_figs_dir,
        ref,
        plum_harm,
        year_list_harm,
        cfg.runlist_legend,
    )
    if not cfg.combine_crops:
        _plot_mgmt_timeseries(
            harms_figs_dir,
            year_list_harm,
            plum_harm_nfert,
            cfg.runlist_legend,
            "nfert",
            plum_harm=plum_harm,
            is_crop=ref["isCrop"],
        )
        _plot_mgmt_timeseries(
            harms_figs_dir,
            year_list_harm,
            plum_harm_irrig,
            cfg.runlist_legend,
            "irrig",
            plum_harm=plum_harm,
            is_crop=ref["isCrop"],
        )

    print("Data arrays prepared for plotting.")
    return {
        "plum_orig_xvyr": plum_orig,
        "plum_harm_xvyr": plum_harm,
        "plum_orig_nfert_xvyr": plum_orig_nfert,
        "plum_orig_irrig_xvyr": plum_orig_irrig,
        "plum_harm_nfert_xvyr": plum_harm_nfert,
        "plum_harm_irrig_xvyr": plum_harm_irrig,
        "harms_figs_dir": harms_figs_dir,
        "year_list_orig": np.array(year_list_orig),
        "year_list_harm": np.array(year_list_harm),
    }


def _to_plumharm_config(cfg: PlumHarmFigsConfig):
    from .plumharm_options import PlumHarmConfig

    return PlumHarmConfig(
        this_dir=cfg.this_dir,
        geodata_dir=cfg.geodata_dir,
        plum_dirs=cfg.plum_dirs,
        harm_dirs=cfg.harm_dirs,
        remap_lu_file=cfg.remap_lu_file,
        remap_cropf_file=cfg.remap_cropf_file,
        remap_nfert_file=cfg.remap_nfert_file,
        base_year=cfg.base_year,
        year1=cfg.year1,
        yearN=cfg.yearN,
        allow_unveg=cfg.allow_unveg,
        conserv_tol_pct=0.2,
        conserv_tol_area=1e3,
        norm2extra=cfg.norm2extra,
        use_latest_plum_mgmt=True,
        inpaint_method=4,
        fix_tiny_negs_tol_m2=1.0,
        save_halfdeg_mat=True,
        save_2deg_mat=True,
        save_halfdeg_txt=False,
        save_2deg_txt=False,
        out_prec=6,
        out_width=1,
        delimiter=" ",
        overwrite=True,
        fancy=False,
        debug_areas=False,
        debug_nfert=False,
        debug_irrig=False,
        debugIJ_2deg=None,
        dbCrop="",
        verbose=True,
        combine_crops=cfg.combine_crops,
        fruitveg_sugar_2oil=cfg.fruitveg_sugar_2oil,
    )


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


def _make_crops_timeseries_fig(
    out_dir: str,
    var_names: Sequence[str],
    year_list_orig: Sequence[int],
    year_list_harm: Sequence[int],
    ts_orig_vyr: np.ndarray,
    ts_harm_vyr: np.ndarray,
    runlist_legend: Sequence[str] | None,
    suffix: str = "area",
    title_word: str = "Area",
    units: str = "m²",
) -> None:
    """Single figure with ceil(N/2)x2 subplots matching MATLAB make_crops_timeseries_fig."""
    if ts_orig_vyr is None or ts_harm_vyr is None:
        return
    n_vars = len(var_names)
    n_rows = int(np.ceil(n_vars / 2))
    n_cols = 2
    n_runs = ts_orig_vyr.shape[-1]
    if runlist_legend is None:
        runlist = [f"run{r+1}" for r in range(n_runs)]
    else:
        runlist = list(runlist_legend)
        if len(runlist) != n_runs:
            runlist = [f"run{r+1}" for r in range(n_runs)]

    fig_h = max(4, 3 * n_rows)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, fig_h), squeeze=False)

    for v in range(n_vars):
        row_idx = v // n_cols
        col_idx = v % n_cols
        ax = axes[row_idx, col_idx]
        for r in range(n_runs):
            ax.plot(year_list_orig, ts_orig_vyr[v, :, r], "--", linewidth=1,
                    label=f"{runlist[r]} orig")
            ax.plot(year_list_harm, ts_harm_vyr[v, :, r], "-", linewidth=1,
                    label=f"{runlist[r]} harm")
        ax.set_title(f"{title_word}: {var_names[v]}")
        ax.set_ylabel(units)
        ax.legend(loc="best", fontsize=7)

    for v in range(n_vars, n_rows * n_cols):
        row_idx = v // n_cols
        col_idx = v % n_cols
        axes[row_idx, col_idx].set_visible(False)

    plt.tight_layout()
    out_path = os.path.join(out_dir, f"timeSeries_{suffix}.pdf")
    plt.savefig(out_path)
    plt.close()


def _plot_delta_maps(
    out_dir: str,
    lu_names: Sequence[str],
    mask_yx: np.ndarray,
    list2map: np.ndarray,
    ny: int,
    nx: int,
    plum_orig: np.ndarray,
    plum_harm: np.ndarray,
    year_list_orig: Sequence[int],
    year_list_harm: Sequence[int],
) -> None:
    y1 = year_list_harm[0]
    yN = year_list_harm[-1]
    i_y1_orig = list(year_list_orig).index(y1)
    i_yN_orig = list(year_list_orig).index(yN)
    i_y1_harm = list(year_list_harm).index(y1)
    i_yN_harm = list(year_list_harm).index(yN)

    for v, lu in enumerate(lu_names):
        if lu not in {"CROPLAND", "PASTURE", "NATURAL"}:
            continue
        orig_delta = plum_orig[:, v, i_yN_orig, 0] - plum_orig[:, v, i_y1_orig, 0]
        harm_delta = plum_harm[:, v, i_yN_harm, 0] - plum_harm[:, v, i_y1_harm, 0]
        orig_map = _vector_to_map(orig_delta, list2map, ny, nx, mask_yx)
        harm_map = _vector_to_map(harm_delta, list2map, ny, nx, mask_yx)

        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        axes[0].imshow(orig_map, cmap="coolwarm")
        axes[0].set_title(f"Orig {lu} delta {y1}-{yN}")
        axes[1].imshow(harm_map, cmap="coolwarm")
        axes[1].set_title(f"Harm {lu} delta {y1}-{yN}")
        for ax in axes:
            ax.axis("off")
        plt.tight_layout()
        out_path = os.path.join(out_dir, f"maps_deltas_{lu}_{y1}-{yN}_beforeAfter.png")
        plt.savefig(out_path, dpi=150)
        plt.close()


def _plot_scatter_baseline_vs_orig(
    out_dir: str,
    ref: Dict[str, np.ndarray],
    plum_orig: np.ndarray,
    plum_harm: np.ndarray,
    year_list_orig: Sequence[int],
    year_list_harm: Sequence[int],
    is2deg: bool,
    ny: int,
    nx: int,
    runlist_legend: Sequence[str] | None = None,
) -> None:
    """Scatter: harm crop fraction vs orig crop fraction at year1 (Hurtt Fig 4)."""
    is_crop = ref.get("isCrop")
    if is_crop is None or not np.any(is_crop):
        return
    y1 = year_list_harm[0]
    i_y1_orig = list(year_list_orig).index(y1)
    i_y1_harm = list(year_list_harm).index(y1)

    land_area = ref["landArea_2deg_YX"] if is2deg else ref["landArea_YX"]
    list2map = ref["list2map_2deg"] if is2deg else ref["list2map"]
    land_area_x = land_area.ravel(order="F")[list2map - 1]

    n_runs = plum_orig.shape[3]
    runlist = list(runlist_legend) if runlist_legend else [f"run{r+1}" for r in range(n_runs)]
    if len(runlist) != n_runs:
        runlist = [f"run{r+1}" for r in range(n_runs)]

    fig, axes = plt.subplots(1, n_runs, figsize=(5 * n_runs, 5), squeeze=False)
    for r in range(n_runs):
        ax = axes[0, r]
        harm_crop = np.sum(plum_harm[:, is_crop, i_y1_harm, r], axis=1) / land_area_x
        orig_crop = np.sum(plum_orig[:, is_crop, i_y1_orig, r], axis=1) / land_area_x
        ax.plot(harm_crop, orig_crop, ".k", markersize=1)
        ax.plot([0, 1], [0, 1], "k--")
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.set_aspect("equal")
        ax.set_xlabel(f"Fraction of gridcell {y1} (Baseline LU)")
        ax.set_ylabel(f"Fraction of gridcell {y1} (PLUM output)")
        if len(runlist) > r:
            ax.set_title(runlist[r])
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "scatter_hurtt2011_fig4.png"), dpi=300)
    plt.close()


def _plot_scatter_orig_vs_harm_delta(
    out_dir: str,
    ref: Dict[str, np.ndarray],
    plum_orig: np.ndarray,
    plum_harm: np.ndarray,
    year_list_orig: Sequence[int],
    year_list_harm: Sequence[int],
    runlist_legend: Sequence[str] | None = None,
) -> None:
    """Scatter: orig vs harm delta for crop and pasture, all runs (Hurtt Fig 5)."""
    is_crop = ref.get("isCrop")
    if is_crop is None or "PASTURE" not in ref["lu_names"]:
        return
    i_past = ref["lu_names"].index("PASTURE")

    y1 = year_list_harm[0]
    yN = year_list_harm[-1]
    i_y1_orig = list(year_list_orig).index(y1)
    i_yN_orig = list(year_list_orig).index(yN)
    i_y1_harm = list(year_list_harm).index(y1)
    i_yN_harm = list(year_list_harm).index(yN)

    land_area_yx = ref["landArea_YX"]
    list2map = ref["list2map"]
    land_area_x = land_area_yx.ravel(order="F")[list2map - 1]
    map_shape = land_area_yx.shape

    n_runs = plum_orig.shape[3]
    runlist = list(runlist_legend) if runlist_legend else [f"run{r+1}" for r in range(n_runs)]
    if len(runlist) != n_runs:
        runlist = [f"run{r+1}" for r in range(n_runs)]

    diff_crop_orig_xr = (
        np.sum(plum_orig[:, is_crop, i_yN_orig, :] - plum_orig[:, is_crop, i_y1_orig, :], axis=1)
        / land_area_x[:, None]
    )
    diff_past_orig_xr = (
        (plum_orig[:, i_past, i_yN_orig, :] - plum_orig[:, i_past, i_y1_orig, :])
        / land_area_x[:, None]
    )
    diff_crop_harm_xr = (
        np.sum(plum_harm[:, is_crop, i_yN_harm, :] - plum_harm[:, is_crop, i_y1_harm, :], axis=1)
        / land_area_x[:, None]
    )
    diff_past_harm_xr = (
        (plum_harm[:, i_past, i_yN_harm, :] - plum_harm[:, i_past, i_y1_harm, :])
        / land_area_x[:, None]
    )

    diff_orig_xrL = np.stack([diff_crop_orig_xr, diff_past_orig_xr], axis=-1)
    diff_harm_xrL = np.stack([diff_crop_harm_xr, diff_past_harm_xr], axis=-1)

    list2map_2deg = ref.get("list2map_2deg")
    diff2_orig: np.ndarray | None = None
    diff2_harm: np.ndarray | None = None
    if list2map_2deg is not None:
        n2 = len(list2map_2deg)
        diff2_orig = np.full((n2, n_runs, 2), np.nan)
        diff2_harm = np.full((n2, n_runs, 2), np.nan)
        for r in range(n_runs):
            for L in range(2):
                diff2_orig[:, r, L] = _aggregate_to_2deg(
                    diff_orig_xrL[:, r, L] * land_area_x, list2map, map_shape,
                    land_area_yx, list2map_2deg,
                )
                diff2_harm[:, r, L] = _aggregate_to_2deg(
                    diff_harm_xrL[:, r, L] * land_area_x, list2map, map_shape,
                    land_area_yx, list2map_2deg,
                )

    gray = (0.65, 0.65, 0.65)
    fig, axes = plt.subplots(2, n_runs, figsize=(5 * n_runs, 10), squeeze=False)
    for r in range(n_runs):
        for L in range(2):
            ax = axes[L, r]
            ax.plot(diff_orig_xrL[:, r, L], diff_harm_xrL[:, r, L], ".",
                    color=gray, markersize=1)
            if diff2_orig is not None:
                ax.plot(diff2_orig[:, r, L], diff2_harm[:, r, L], ".k", markersize=3)
            ax.plot([-1, 1], [-1, 1], "--k")
            ax.set_xlim([-1, 1])
            ax.set_ylim([-1, 1])
            ax.set_aspect("equal")
            ax.set_xlabel("\u0394 gridcell fraction (original)")
            ax.set_ylabel("\u0394 gridcell fraction (harmonized)")
            if L == 0:
                ax.set_title(runlist[r])
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "scatter_hurtt2011_fig5.png"), dpi=300)
    plt.close()


def _vector_to_map(vec: np.ndarray, list2map: np.ndarray, ny: int, nx: int, mask_yx: np.ndarray):
    out = np.full((ny, nx), np.nan)
    for idx, val in zip(list2map, vec):
        i, j = np.unravel_index(idx - 1, (ny, nx), order="F")
        out[i, j] = val
    out[mask_yx] = np.nan
    return out


def _aggregate_to_2deg(
    weighted_x: np.ndarray,
    list2map: np.ndarray,
    map_shape: Tuple[int, int],
    land_area_yx: np.ndarray,
    list2map_2deg: np.ndarray,
) -> np.ndarray:
    """Aggregate half-degree area-weighted values to 2-degree fraction means."""
    full_map = np.zeros(map_shape)
    for idx, val in zip(list2map, weighted_x):
        i, j = np.unravel_index(idx - 1, map_shape, order="F")
        if not np.isnan(val):
            full_map[i, j] = val

    agg = full_map[:, 0::4] + full_map[:, 1::4] + full_map[:, 2::4] + full_map[:, 3::4]
    agg = agg[0::4, :] + agg[1::4, :] + agg[2::4, :] + agg[3::4, :]

    area_agg = (
        land_area_yx[:, 0::4] + land_area_yx[:, 1::4]
        + land_area_yx[:, 2::4] + land_area_yx[:, 3::4]
    )
    area_agg = area_agg[0::4, :] + area_agg[1::4, :] + area_agg[2::4, :] + area_agg[3::4, :]

    with np.errstate(invalid="ignore", divide="ignore"):
        result_yx = np.where(area_agg > 0, agg / area_agg, np.nan)

    return result_yx.ravel(order="F")[list2map_2deg - 1]


def _write_geotiffs(
    out_dir: str,
    ref: Dict[str, Any],
    plum_orig: np.ndarray,
    plum_harm: np.ndarray,
    year_list_orig: Sequence[int],
    year_list_harm: Sequence[int],
    this_ver: str,
    runlist_legend: Sequence[str] | None,
) -> None:
    try:
        import rasterio
        from rasterio.transform import from_origin
    except Exception as exc:
        raise RuntimeError("rasterio is required for GeoTIFF export.") from exc

    is2deg = this_ver == "2deg."
    ny, nx = (90, 180) if is2deg else (360, 720)
    list2map = ref["list2map_2deg"] if is2deg else ref["list2map"]
    mask_yx = ref["mask_YX"]
    res = 2.0 if is2deg else 0.5
    transform = from_origin(-180, 90, res, res)

    y1 = year_list_harm[0]
    yN = year_list_harm[-1]
    i_y1_orig = list(year_list_orig).index(y1)
    i_yN_orig = list(year_list_orig).index(yN)
    i_y1_harm = list(year_list_harm).index(y1)
    i_yN_harm = list(year_list_harm).index(yN)

    n_runs = plum_orig.shape[3]
    if runlist_legend is None:
        runlist = [f"run{r+1}" for r in range(n_runs)]
    else:
        runlist = list(runlist_legend)

    lu_short = {"CROPLAND": "crop", "PASTURE": "past", "NATURAL": "ntrl"}

    for lu in ["CROPLAND", "PASTURE", "NATURAL"]:
        if lu not in ref["lu_names"]:
            continue
        v = ref["lu_names"].index(lu)
        for r in range(n_runs):
            orig_delta = plum_orig[:, v, i_yN_orig, r] - plum_orig[:, v, i_y1_orig, r]
            harm_delta = plum_harm[:, v, i_yN_harm, r] - plum_harm[:, v, i_y1_harm, r]
            orig_map = _vector_to_map(orig_delta, list2map, ny, nx, mask_yx)
            harm_map = _vector_to_map(harm_delta, list2map, ny, nx, mask_yx)

            for tag, data in [("orig", orig_map), ("harm", harm_map)]:
                out_path = os.path.join(
                    out_dir, f"D{lu_short[lu]}_{y1}-{yN}_{runlist[r]}_{tag}.tif"
                )
                with rasterio.open(
                    out_path,
                    "w",
                    driver="GTiff",
                    height=data.shape[0],
                    width=data.shape[1],
                    count=1,
                    dtype=data.dtype,
                    crs="EPSG:4326",
                    transform=transform,
                    nodata=-999,
                ) as dst:
                    dst.write(np.flipud(data).astype(data.dtype), 1)


def _harm_by_numbers(
    out_dir: str,
    ref: Dict[str, Any],
    plum_orig: np.ndarray,
    plum_harm: np.ndarray,
    year_list_orig: Sequence[int],
    year_list_harm: Sequence[int],
    runlist_legend: Sequence[str] | None,
    combine_crops: bool,
) -> None:
    incl_years = [y for y in year_list_orig if 2010 <= y <= 2100]
    yi_orig = [list(year_list_orig).index(y) for y in incl_years]
    yi_harm = [list(year_list_harm).index(y) for y in incl_years if y in year_list_harm]
    if len(yi_orig) != len(incl_years) or len(yi_harm) != len(incl_years):
        raise RuntimeError("Mismatch in included years for harm_by_numbers.")

    y1 = min(incl_years)
    yN = max(incl_years)
    tmp_lu_list = ["NATURAL", "CROPLAND", "PASTURE"]

    biome_id_x, countries_x, countries_key = _load_region_masks(ref)
    lons, lats = _lonlat_from_list2map(ref["list2map"], ref["landArea_YX"].shape)
    land_area_x = ref["landArea_YX"].ravel(order="F")[ref["list2map"] - 1]

    harm_focus_regions = _harm_focus_regions(biome_id_x, countries_x, countries_key)
    list_regions = [row[4] for row in harm_focus_regions]
    list_super_regs = list(dict.fromkeys([row[1] for row in harm_focus_regions]))

    template = os.path.join(ref.get("template_dir", ""), "harm_by_numbers.template.xlsx")
    template_rel = os.path.join(ref.get("template_dir", ""), "harm_by_numbers.template.relLandArea.xlsx")
    if not os.path.exists(template):
        template = None
    if not os.path.exists(template_rel):
        template_rel = None

    n_runs = plum_orig.shape[3]
    if runlist_legend is None:
        runlist = [f"run{r+1}" for r in range(n_runs)]
    else:
        runlist = list(runlist_legend)

    for lu in tmp_lu_list:
        outfile = os.path.join(out_dir, f"harm_by_numbers.{y1}-{yN}.{lu}.xlsx")
        if combine_crops:
            outfile = outfile.replace(lu, f"combCrops.{lu}")
        outfile_rel = outfile.replace(".xlsx", ".relLandArea.xlsx")
        if template:
            shutil.copyfile(template, outfile)
        if template_rel:
            shutil.copyfile(template_rel, outfile_rel)

        if lu == "CROPLAND":
            v_mask = ref["isCrop"]
        else:
            v_mask = np.array([v == lu for v in ref["lu_names"]])

        for r in range(n_runs):
            land_area_by_reg = None
            table_orig_out = []
            table_harm_out = []
            table_orig_rel_out = []
            table_harm_rel_out = []
            land_area_by_reg_out = []

            for s_idx, _sreg in enumerate(list_super_regs):
                table_orig, table_orig_rel, land_area_by_reg = _iterate_superreg(
                    np.sum(plum_orig[:, v_mask, :, r][:, :, yi_orig], axis=1),
                    s_idx,
                    "orig",
                    list_super_regs,
                    list_regions,
                    harm_focus_regions,
                    biome_id_x,
                    countries_x,
                    countries_key,
                    lats,
                    lons,
                    ref["list2map"],
                    land_area_x,
                )
                table_harm, table_harm_rel, _ = _iterate_superreg(
                    np.sum(plum_harm[:, v_mask, :, r][:, :, yi_harm], axis=1),
                    s_idx,
                    "harm",
                    list_super_regs,
                    list_regions,
                    harm_focus_regions,
                    biome_id_x,
                    countries_x,
                    countries_key,
                    lats,
                    lons,
                    ref["list2map"],
                    land_area_x,
                )
                table_orig_out.append(table_orig)
                table_orig_rel_out.append(table_orig_rel)
                table_harm_out.append(table_harm)
                table_harm_rel_out.append(table_harm_rel)
                land_area_by_reg_out.append(land_area_by_reg)

            table_orig_out = pd.concat(table_orig_out, ignore_index=True)
            table_harm_out = pd.concat(table_harm_out, ignore_index=True)
            table_orig_rel_out = pd.concat(table_orig_rel_out, ignore_index=True)
            table_harm_rel_out = pd.concat(table_harm_rel_out, ignore_index=True)
            land_area_by_reg_out = np.concatenate(land_area_by_reg_out)

            table_orig_rel_land = table_orig_out.copy()
            table_harm_rel_land = table_harm_out.copy()
            for col in ["Gain", "Loss", "Net"]:
                table_orig_rel_land[col] = table_orig_out[col] / land_area_by_reg_out
                table_harm_rel_land[col] = table_harm_out[col] / land_area_by_reg_out

            sheet_prefix = runlist[r].replace(".", "_")
            _write_excel_sheets(outfile, sheet_prefix, table_orig_out, table_harm_out, table_orig_rel_out, table_harm_rel_out)
            _write_excel_sheets(outfile_rel, sheet_prefix, table_orig_rel_land, table_harm_rel_land, table_orig_rel_out, table_harm_rel_out)


def _write_excel_sheets(
    filename: str,
    sheet_prefix: str,
    table_orig: pd.DataFrame,
    table_harm: pd.DataFrame,
    table_orig_rel: pd.DataFrame,
    table_harm_rel: pd.DataFrame,
) -> None:
    mode = "a" if os.path.exists(filename) else "w"
    with pd.ExcelWriter(filename, engine="openpyxl", mode=mode, if_sheet_exists="replace") as writer:
        table_orig.to_excel(writer, sheet_name=f"{sheet_prefix}_o2", index=False)
        table_harm.to_excel(writer, sheet_name=f"{sheet_prefix}_h2", index=False)
        table_orig_rel.to_excel(writer, sheet_name=f"{sheet_prefix}_o2Ry1", index=False)
        table_harm_rel.to_excel(writer, sheet_name=f"{sheet_prefix}_h2Ry1", index=False)


def _iterate_superreg(
    area_x1y: np.ndarray,
    s_idx: int,
    _orig_or_harm: str,
    list_super_regs: Sequence[str],
    list_regions: Sequence[str],
    harm_focus_regions: Sequence[Sequence[Any]],
    biome_id_x: np.ndarray,
    countries_x: np.ndarray,
    countries_key: pd.DataFrame,
    lats: np.ndarray,
    lons: np.ndarray,
    list2map: np.ndarray,
    land_area_x: np.ndarray,
) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    this_super = list_super_regs[s_idx]
    region_in_super = [row[1] == this_super for row in harm_focus_regions]
    subset = [row for row, ok in zip(harm_focus_regions, region_in_super) if ok]

    list_gain: List[float] = []
    list_loss: List[float] = []
    list_net: List[float] = []
    list_land: List[float] = []
    list_orig: List[float] = []
    list_gain_rel: List[float] = []
    list_loss_rel: List[float] = []
    list_net_rel: List[float] = []
    list_subreg: List[str] = []
    list_biome: List[str] = []
    list_super_out: List[str] = []
    super_gain = super_loss = super_net = super_orig = 0.0

    biome_total_indices: Dict[str, int] = {}
    biome_orig_accum: Dict[str, float] = {}

    for s2, row in enumerate(subset):
        incl_regions = row[0]
        this_biome = row[2]
        tmp_incl_countries = row[3]
        this_reg_name = row[4]

        do_biome_total = sum(1 for r2 in subset if r2[2] == this_biome) > 1
        if do_biome_total and this_biome not in biome_total_indices:
            bt_idx = len(list_gain)
            biome_total_indices[this_biome] = bt_idx
            biome_orig_accum[this_biome] = 0.0
            list_gain.append(0.0)
            list_loss.append(0.0)
            list_net.append(0.0)
            list_land.append(0.0)
            list_orig.append(0.0)
            list_gain_rel.append(0.0)
            list_loss_rel.append(0.0)
            list_net_rel.append(0.0)
            list_subreg.append("TOTAL (BIOME)")
            list_biome.append(this_biome)
            list_super_out.append(this_super)

        is_cell_in_reg = np.zeros_like(biome_id_x, dtype=bool)
        if isinstance(tmp_incl_countries, list) and len(tmp_incl_countries) > 0:
            incl_codes = countries_key.loc[
                countries_key["Country"].isin(tmp_incl_countries), "numCode"
            ].to_numpy()
            is_cell_in_reg |= np.isin(biome_id_x, incl_regions) & np.isin(countries_x, incl_codes)
        elif isinstance(tmp_incl_countries, str) and "&" in tmp_incl_countries:
            bound, news = tmp_incl_countries.split("&")
            bound = float(bound)
            if news.lower() == "n":
                is_cell_in_reg |= np.isin(biome_id_x, incl_regions) & (lats >= bound)
            elif news.lower() == "s":
                is_cell_in_reg |= np.isin(biome_id_x, incl_regions) & (lats <= bound)
            elif news.lower() == "e":
                is_cell_in_reg |= np.isin(biome_id_x, incl_regions) & (lons >= bound)
            elif news.lower() == "w":
                is_cell_in_reg |= np.isin(biome_id_x, incl_regions) & (lons <= bound)
            else:
                raise RuntimeError("Unrecognized geo restriction")
        else:
            is_cell_in_reg |= np.isin(biome_id_x, incl_regions)

        gain, loss, orig = _get_region_numbers(area_x1y[is_cell_in_reg, :])
        net = gain + loss
        super_gain += gain
        super_loss += loss
        super_net += net
        super_orig += orig
        land_area = np.sum(land_area_x[is_cell_in_reg])

        list_gain.append(gain)
        list_loss.append(loss)
        list_net.append(net)
        list_land.append(land_area)
        list_orig.append(orig)
        list_gain_rel.append(gain / orig if orig != 0 else np.nan)
        list_loss_rel.append(loss / orig if orig != 0 else np.nan)
        list_net_rel.append(net / orig if orig != 0 else np.nan)
        list_subreg.append(this_reg_name)
        list_biome.append(this_biome)
        list_super_out.append(this_super)

        if do_biome_total:
            bt_idx = biome_total_indices[this_biome]
            list_gain[bt_idx] += gain
            list_loss[bt_idx] += loss
            list_net[bt_idx] += net
            list_land[bt_idx] += land_area
            biome_orig_accum[this_biome] += orig
            is_last = s2 == max(
                j for j, r2 in enumerate(subset) if r2[2] == this_biome
            )
            if is_last:
                b_orig = biome_orig_accum[this_biome]
                list_gain_rel[bt_idx] = list_gain[bt_idx] / b_orig if b_orig != 0 else np.nan
                list_loss_rel[bt_idx] = list_loss[bt_idx] / b_orig if b_orig != 0 else np.nan
                list_net_rel[bt_idx] = list_net[bt_idx] / b_orig if b_orig != 0 else np.nan

    if len(subset) > 1:
        list_gain.insert(0, np.nansum(list_gain))
        list_loss.insert(0, np.nansum(list_loss))
        list_net.insert(0, np.nansum(list_net))
        list_land.insert(0, np.nansum(list_land))
        list_orig.insert(0, 0.0)
        list_gain_rel.insert(0, super_gain / super_orig if super_orig != 0 else np.nan)
        list_loss_rel.insert(0, super_loss / super_orig if super_orig != 0 else np.nan)
        list_net_rel.insert(0, super_net / super_orig if super_orig != 0 else np.nan)
        list_subreg.insert(0, "TOTAL (SUPREG)")
        list_biome.insert(0, "TOTAL")
        list_super_out.insert(0, this_super)

    table_out = pd.DataFrame(
        {
            "Super-region": list_super_out,
            "Biome": list_biome,
            "Sub-region": list_subreg,
            "Gain": list_gain,
            "Loss": list_loss,
            "Net": list_net,
        }
    )
    table_out_rel = pd.DataFrame(
        {
            "Super-region": list_super_out,
            "Biome": list_biome,
            "Sub-region": list_subreg,
            "Gain": list_gain_rel,
            "Loss": list_loss_rel,
            "Net": list_net_rel,
        }
    )
    return table_out, table_out_rel, np.array(list_land)


def _get_region_numbers(area_x1y: np.ndarray) -> Tuple[float, float, float]:
    orig_x = area_x1y[:, 0]
    diff_x = area_x1y[:, -1] - orig_x
    gain = np.sum(diff_x[diff_x > 0])
    loss = np.sum(diff_x[diff_x < 0])
    orig = np.sum(orig_x)
    return gain, loss, orig


def _load_region_masks(ref: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    try:
        import rasterio
    except Exception as exc:
        raise RuntimeError("rasterio is required to read biome GeoTIFF.") from exc

    biome_path = os.path.join(ref["geodata_dir"], "wwf_terr_ecos_dissolveBiome_halfDeg_id.tif")
    with rasterio.open(biome_path) as src:
        biome_id_yx = src.read(1)
    biome_id_yx = np.flipud(biome_id_yx)
    biome_id_yx = np.where(biome_id_yx < 0, np.nan, biome_id_yx)

    countries_path = os.path.join(
        ref["geodata_dir"], "country_boundaries", "country_boundaries62892.noNeg99.extrapd.asc"
    )
    countries_yx = np.loadtxt(countries_path, skiprows=6)
    countries_yx = np.flipud(countries_yx)
    countries_yx = np.where(countries_yx <= 0, np.nan, countries_yx)

    countries_key = pd.read_csv(
        os.path.join(ref["geodata_dir"], "country_boundaries", "country_boundaries_codes4.csv")
    )
    list2map = ref["list2map"]
    biome_id_x = biome_id_yx.ravel(order="F")[list2map - 1]
    countries_x = countries_yx.ravel(order="F")[list2map - 1]
    return biome_id_x, countries_x, countries_key


def _lonlat_from_list2map(list2map: np.ndarray, map_shape: Tuple[int, int]) -> Tuple[np.ndarray, np.ndarray]:
    nlat, nlon = map_shape
    xres = 360 / nlon
    yres = 180 / nlat
    lons = np.arange(-180 + xres / 2, 180 - xres / 2 + 1e-9, xres)
    lats = np.arange(-90 + yres / 2, 90 - yres / 2 + 1e-9, yres)
    lons_map = np.tile(lons, (len(lats), 1))
    lats_map = np.tile(lats.reshape(-1, 1), (1, len(lons)))
    lons_x = lons_map.ravel(order="F")[list2map - 1]
    lats_x = lats_map.ravel(order="F")[list2map - 1]
    return lons_x, lats_x


def _harm_focus_regions(
    biome_id_x: np.ndarray, countries_x: np.ndarray, countries_key: pd.DataFrame
) -> List[List[Any]]:
    china_code = countries_key.loc[countries_key["Country"] == "China", "numCode"].iloc[0]
    china_biomes = set(np.unique(biome_id_x[countries_x == china_code]))
    china_other = list(china_biomes - {480, 487, 342, 348, 662})
    world_biomes = [v for v in np.unique(biome_id_x) if not np.isnan(v)]
    return [
        [world_biomes, "World", "World", [], "World"],
        [1095, "Amazon", "Trop. rainforest", [], "Amazon"],
        [429, "N. America", "Temp. grassland", [], "Great Plains"],
        [437, "N. America", "Temp. forest", [], "E US mixed for"],
        [306, "N. America", "Temp. forest", [], "U. Midw US br/mix for"],
        [477, "N. America", "Temp. forest", [], "E US conif for"],
        [436, "N. America", "Temp. forest", [], "Texarkana conif for"],
        [229, "Alaska", "Tundra", [], "Alaskan forest"],
        [[62, 63, 68, 72, 78, 79, 80, 87, 95, 96, 97, 101, 102, 109, 110, 111, 112, 113,
          124, 125, 126, 127, 131, 133, 134, 135, 136, 139, 143, 144, 150, 163, 175, 202,
          213, 220, 229], "Alaska", "Bor. forest", ["United States of America"], "Alaskan tundra"],
        [950, "Sub-Sah. Afr.", "Trop. rainforest", [], "C Afr rainfor."],
        [1251, "Sub-Sah. Afr.", "Savanna", "0&N", "N Afr savanna"],
        [1251, "Sub-Sah. Afr.", "Savanna", "0&S", "S Afr savanna"],
        [811, "South Asia", "Desert/xeric", ["India", "Pakistan", "Afghanistan"], "S Asia xeric/desert"],
        [661, "South Asia", "(Sub)trop. dry for.", [], "C Ind subt dry for"],
        [712, "South Asia", "(Sub)trop. dry for.", [], "S Ind subt dry for"],
        [766, "South Asia", "(Sub)trop. dry for.", [], "S Ind scrub for"],
        [808, "South Asia", "(Sub)trop. dry for.", [], "SriL subt dry for"],
        [708, "South Asia", "(Sub)trop. wet for.", [], "W Ind subt wet for"],
        [830, "South Asia", "(Sub)trop. wet for.", [], "SriL subt wet for"],
        [597, "South Asia", "(Sub)trop. wet for.", [], "C Ind subt wet for"],
        [[551, 567, 632, 638], "South Asia", "(Sub)trop. wet for.", [], "E Ind subt wet for"],
        [662, "South Asia", "(Sub)trop. wet for.", ["India", "Bangladesh"], "NWInd+Bangl subt wet for"],
        [480, "East Asia", "Temp. forest", [], "E Asia temp for"],
        [487, "East Asia", "Montane gr/shr", [], "Tibetan Plat. steppe"],
        [342, "East Asia", "Desert/xeric", [], "E Asia xeric/desert"],
        [348, "East Asia", "Temp. grass/sav/shr", [], "E Asia temp grass"],
        [480, "China", "Temp. forest", ["China"], "E China temp for"],
        [487, "China", "Montane gr/shr", ["China"], "China Tib. Plat. steppe"],
        [342, "China", "Desert/xeric", ["China"], "China xeric/desert"],
        [348, "China", "Temp. grass/sav/shr", ["China"], "China temp grass"],
        [662, "China", "(Sub)trop. wet for.", ["China"], "China subt wet for"],
        [china_other, "China", "China other", ["China"], "China other"],
        [[132, 153, 169, 171, 252, 313], "Europe+Nafr", "Temp. forest", [], "Eur temp br/mix for"],
        [[168, 287, 288, 347], "Europe+Nafr", "Temp. forest", [], "Eur temp conif for"],
        [349, "Europe+Nafr", "Mediterranean", [], "Mediterr. mediterr."],
    ]


def _plot_harm_effect_nonagri(
    out_dir: str,
    ref: Dict[str, Any],
    plum_orig: np.ndarray,
    plum_harm: np.ndarray,
    year_list_orig: Sequence[int],
    runlist_legend: Sequence[str] | None,
    legend_loc: str,
) -> None:
    if "NATURAL" not in ref["lu_names"]:
        return
    idx = ref["lu_names"].index("NATURAL")
    orig_incr = plum_orig[:, idx, 1:, :] - plum_orig[:, idx, :-1, :]
    harm_incr = plum_harm[:, idx, 1:, :] - plum_harm[:, idx, :-1, :]
    orig_incr[orig_incr < 0] = 0
    harm_incr[harm_incr < 0] = 0
    orig_sum = np.sum(orig_incr, axis=0)
    harm_sum = np.sum(harm_incr, axis=0)
    harm_effect = harm_sum - orig_sum
    x = np.array(year_list_orig[1:])

    plt.figure(figsize=(12, 3))
    plt.plot(x, harm_effect, linewidth=1)
    for r in range(harm_effect.shape[1]):
        coeff = np.polyfit(x, harm_effect[:, r], 1)
        plt.plot(x, np.polyval(coeff, x), "--", linewidth=2)
    if runlist_legend is not None:
        plt.legend(runlist_legend, loc=legend_loc)
    plt.title("Time series of harmonization effect on change in non-agri area")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "timeSeries_harm_effect_on_change_in_nonagri_area.pdf"))
    plt.close()


def _plot_bundle_timeseries(
    out_dir: str,
    ref: Dict[str, Any],
    plum_harm: np.ndarray,
    year_list_harm: Sequence[int],
    runlist_legend: Sequence[str] | None,
) -> None:
    n_runs = plum_harm.shape[3]
    runlist = list(runlist_legend) if runlist_legend else [f"run{r+1}" for r in range(n_runs)]

    # Land use
    plt.figure(figsize=(10, 4))
    for v, name in enumerate(ref["lu_names"]):
        series = np.nansum(plum_harm[:, v, :, 0], axis=0)
        plt.plot(year_list_harm, series, label=name)
    plt.title("Time series land use")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "timeSeries_landUse.pdf"))
    plt.close()

    # Crops (from isCrop mask in LU list)
    if "isCrop" in ref:
        plt.figure(figsize=(10, 4))
        for v, name in enumerate([n for n, is_c in zip(ref["lu_names"], ref["isCrop"]) if is_c]):
            idx = ref["lu_names"].index(name)
            series = np.nansum(plum_harm[:, idx, :, 0], axis=0)
            plt.plot(year_list_harm, series, label=name)
        plt.title("Time series crops")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, "timeSeries_crops.pdf"))
        plt.close()

    # Placeholders for nfert/irrig handled by _make_crops_timeseries_fig
    if os.path.exists(os.path.join(out_dir, "timeSeries_nfert.pdf")):
        pass
    if os.path.exists(os.path.join(out_dir, "timeSeries_irrig.pdf")):
        pass


def _plot_mgmt_timeseries(
    out_dir: str,
    year_list_harm: Sequence[int],
    harm_xvyr: np.ndarray | None,
    runlist_legend: Sequence[str] | None,
    suffix: str,
    plum_harm: np.ndarray | None = None,
    is_crop: np.ndarray | None = None,
) -> None:
    if harm_xvyr is None:
        return
    n_runs = harm_xvyr.shape[3]
    runlist = list(runlist_legend) if runlist_legend else [f"run{r+1}" for r in range(n_runs)]
    plt.figure(figsize=(10, 4))
    for r in range(n_runs):
        if plum_harm is not None and is_crop is not None:
            crop_area = plum_harm[:, is_crop, :, r]
            weighted = harm_xvyr[:, :, :, r] * crop_area
            series = np.nansum(weighted, axis=(0, 1))
            if suffix == "nfert":
                series = series * 1e-3 * 1e-6
        else:
            series = np.nansum(harm_xvyr[:, :, :, r], axis=(0, 1))
        plt.plot(year_list_harm, series, label=runlist[r])
    plt.title(f"Time series {suffix}")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"timeSeries_{suffix}.pdf"))
    plt.close()


def _plot_three_year_maps(
    out_dir: str,
    ref: Dict[str, Any],
    plum_orig: np.ndarray,
    plum_harm: np.ndarray,
    year_list_orig: Sequence[int],
    year_list_harm: Sequence[int],
    three_years: Sequence[int],
    runlist_legend: Sequence[str] | None,
    is2deg: bool,
    ny: int,
    nx: int,
) -> None:
    list2map = ref["list2map_2deg"] if is2deg else ref["list2map"]
    mask_yx = ref["mask_YX"]
    runlist = list(runlist_legend) if runlist_legend else ["run1"]

    for r in range(plum_orig.shape[3]):
        run_name = runlist[r] if r < len(runlist) else f"run{r+1}"
        for lu in ["CROPLAND", "PASTURE", "NATURAL"]:
            if lu not in ref["lu_names"]:
                continue
            v = ref["lu_names"].index(lu)
            fig, axes = plt.subplots(2, len(three_years), figsize=(4 * len(three_years), 6))
            for j, y in enumerate(three_years):
                i_orig = list(year_list_orig).index(y)
                i_harm = list(year_list_harm).index(y)
                orig_map = _vector_to_map(plum_orig[:, v, i_orig, r], list2map, ny, nx, mask_yx)
                harm_map = _vector_to_map(plum_harm[:, v, i_harm, r], list2map, ny, nx, mask_yx)
                axes[0, j].imshow(orig_map, cmap="viridis")
                axes[0, j].set_title(f"Orig {y}")
                axes[1, j].imshow(harm_map, cmap="viridis")
                axes[1, j].set_title(f"Harm {y}")
                axes[0, j].axis("off")
                axes[1, j].axis("off")
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, f"maps_{lu}_{three_years}_{run_name}.png"))
            plt.close()


def _plot_change_maps(
    out_dir: str,
    ref: Dict[str, Any],
    plum_orig: np.ndarray,
    plum_harm: np.ndarray,
    year_list_orig: Sequence[int],
    year_list_harm: Sequence[int],
    three_years: Sequence[int],
    runlist_legend: Sequence[str] | None,
    is2deg: bool,
    ny: int,
    nx: int,
) -> None:
    list2map = ref["list2map_2deg"] if is2deg else ref["list2map"]
    mask_yx = ref["mask_YX"]
    runlist = list(runlist_legend) if runlist_legend else ["run1"]

    if len(three_years) < 2:
        return
    pairs = list(zip(three_years[:-1], three_years[1:]))
    for r in range(plum_orig.shape[3]):
        run_name = runlist[r] if r < len(runlist) else f"run{r+1}"
        for lu in ["CROPLAND", "PASTURE", "NATURAL"]:
            if lu not in ref["lu_names"]:
                continue
            v = ref["lu_names"].index(lu)
            fig, axes = plt.subplots(2, len(pairs), figsize=(4 * len(pairs), 6))
            for j, (y1, y2) in enumerate(pairs):
                i1o = list(year_list_orig).index(y1)
                i2o = list(year_list_orig).index(y2)
                i1h = list(year_list_harm).index(y1)
                i2h = list(year_list_harm).index(y2)
                orig_map = _vector_to_map(
                    plum_orig[:, v, i2o, r] - plum_orig[:, v, i1o, r], list2map, ny, nx, mask_yx
                )
                harm_map = _vector_to_map(
                    plum_harm[:, v, i2h, r] - plum_harm[:, v, i1h, r], list2map, ny, nx, mask_yx
                )
                axes[0, j].imshow(orig_map, cmap="coolwarm")
                axes[0, j].set_title(f"Orig {y1}-{y2}")
                axes[1, j].imshow(harm_map, cmap="coolwarm")
                axes[1, j].set_title(f"Harm {y1}-{y2}")
                axes[0, j].axis("off")
                axes[1, j].axis("off")
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, f"mapsChgs_{lu}_{three_years}_{run_name}.png"))
            plt.close()


def _plot_diff_maps_three_years(
    out_dir: str,
    ref: Dict[str, Any],
    plum_orig: np.ndarray,
    plum_harm: np.ndarray,
    year_list_orig: Sequence[int],
    year_list_harm: Sequence[int],
    three_years: Sequence[int],
    runlist_legend: Sequence[str] | None,
    is2deg: bool,
    ny: int,
    nx: int,
) -> None:
    """Maps of (harm - orig) at each of three_years, per LU and per run."""
    list2map = ref["list2map_2deg"] if is2deg else ref["list2map"]
    mask_yx = ref["mask_YX"]
    runlist = list(runlist_legend) if runlist_legend else ["run1"]

    for r in range(plum_orig.shape[3]):
        run_name = runlist[r] if r < len(runlist) else f"run{r+1}"
        for lu in ["CROPLAND", "PASTURE", "NATURAL"]:
            if lu not in ref["lu_names"]:
                continue
            v = ref["lu_names"].index(lu)
            fig, axes = plt.subplots(1, len(three_years),
                                     figsize=(4 * len(three_years), 3))
            if len(three_years) == 1:
                axes = [axes]
            for j, y in enumerate(three_years):
                if y not in year_list_orig or y not in year_list_harm:
                    continue
                i_orig = list(year_list_orig).index(y)
                i_harm = list(year_list_harm).index(y)
                diff = plum_harm[:, v, i_harm, r] - plum_orig[:, v, i_orig, r]
                diff_map = _vector_to_map(diff, list2map, ny, nx, mask_yx)
                ax = axes[j] if isinstance(axes, list) else axes[j]
                im = ax.imshow(diff_map, cmap="coolwarm")
                vmax = np.nanmax(np.abs(diff_map))
                if vmax > 0:
                    im.set_clim(-vmax, vmax)
                ax.set_title(f"Harm-Orig, {run_name}: {lu}, {y}")
                ax.axis("off")
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            plt.tight_layout()
            yrs_str = "-".join(str(y) for y in three_years)
            plt.savefig(os.path.join(
                out_dir,
                f"maps_diffs_harmMinusOrig_{lu}_{yrs_str}_{run_name}.png",
            ), dpi=150)
            plt.close()


def _plot_diff_maps_one_year_all_runs(
    out_dir: str,
    ref: Dict[str, Any],
    plum_orig: np.ndarray,
    plum_harm: np.ndarray,
    year_list_orig: Sequence[int],
    year_list_harm: Sequence[int],
    year: int,
    runlist_legend: Sequence[str] | None,
    is2deg: bool,
    ny: int,
    nx: int,
) -> None:
    """Maps of (harm - orig) at a single year with one subplot per run."""
    list2map = ref["list2map_2deg"] if is2deg else ref["list2map"]
    mask_yx = ref["mask_YX"]
    n_runs = plum_orig.shape[3]
    runlist = list(runlist_legend) if runlist_legend else [f"run{r+1}" for r in range(n_runs)]
    if len(runlist) != n_runs:
        runlist = [f"run{r+1}" for r in range(n_runs)]

    if year not in year_list_orig or year not in year_list_harm:
        return
    i_orig = list(year_list_orig).index(year)
    i_harm = list(year_list_harm).index(year)

    land_area = ref["landArea_2deg_YX"] if is2deg else ref["landArea_YX"]
    land_area_x = land_area.ravel(order="F")[list2map - 1]

    n_lu = plum_orig.shape[1]
    n_cols = int(np.ceil(np.sqrt(n_runs)))
    n_rows = int(np.ceil(n_runs / n_cols))

    for lu in ["CROPLAND", "PASTURE", "NATURAL"]:
        if lu not in ref["lu_names"]:
            continue
        v = ref["lu_names"].index(lu)
        fig, axes = plt.subplots(n_rows, n_cols,
                                 figsize=(5 * n_cols, 4 * n_rows), squeeze=False)
        global_vmax = 0.0
        maps_list = []
        for r in range(n_runs):
            diff = plum_harm[:, v, i_harm, r] - plum_orig[:, v, i_orig, r]
            diff_frac = 100.0 * diff / land_area_x
            diff_frac[land_area_x == 0] = np.nan
            diff_map = _vector_to_map(diff_frac, list2map, ny, nx, mask_yx)
            maps_list.append(diff_map)
            vmax = np.nanmax(np.abs(diff_map))
            if vmax > global_vmax:
                global_vmax = vmax
        if global_vmax == 0:
            global_vmax = 1.0

        for r in range(n_runs):
            ri, ci = divmod(r, n_cols)
            ax = axes[ri, ci]
            im = ax.imshow(maps_list[r], cmap="coolwarm", vmin=-global_vmax, vmax=global_vmax)
            ax.set_title(runlist[r])
            ax.axis("off")
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        for r in range(n_runs, n_rows * n_cols):
            ri, ci = divmod(r, n_cols)
            axes[ri, ci].set_visible(False)

        fig.suptitle(f"Harm-Orig (%): {lu}, {year}", fontsize=14, fontweight="bold")
        plt.tight_layout()
        plt.savefig(os.path.join(
            out_dir, f"maps_diffs_harmMinusOrig_{lu}_{year}_allRuns.png",
        ), dpi=150)
        plt.close()
