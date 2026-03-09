"""Management (nfert/irrig) harmonization (partial, MATLAB parity)."""
from __future__ import annotations

from typing import Tuple

import numpy as np

from .plumharm_checks import check_bad_vals, check_mgmt_conservation
from .plumharm_debug import debug_out_deltas
from .plumharm_ring_redist import ring_redistribute_mgmt


def harmonize_mgmt(*args, **kwargs):
    """Harmonize management inputs (PLUMharm.m mgmt section)."""
    (
        out_y0_2deg_nfert,
        out_y0_2deg_irrig,
        out_y0_2deg_agri,
        in_y0_2deg_nfert,
        in_y1_2deg_nfert,
        in_y0_2deg_irrig,
        in_y1_2deg_irrig,
        in_y0_2deg_agri,
        in_y1_2deg_agri,
        out_y1_2deg_agri,
        max_nfert_y1,
        lpjg_crops,
        conserv_tol_pct,
        debug_ij,
        db_crop,
        debug_nfert,
        debug_irrig,
        out_prec,
    ) = args

    

    db_crop_idx = None
    if db_crop:
        try:
            db_crop_idx = int(db_crop)
        except (TypeError, ValueError):
            db_crop_idx = lpjg_crops.index(db_crop) if db_crop in lpjg_crops else None

    unm_nfert, mid_nfert = get_unmet_mgmt(
        out_y0_2deg_nfert,
        out_y0_2deg_agri,
        in_y0_2deg_nfert,
        in_y1_2deg_nfert,
        in_y0_2deg_agri,
        in_y1_2deg_agri,
        out_y1_2deg_agri,
        max_nfert_y1,
        debug_ij,
        "Nfert",
        db_crop_idx,
    )
    unm_irrig, mid_irrig = get_unmet_mgmt(
        out_y0_2deg_irrig,
        out_y0_2deg_agri,
        in_y0_2deg_irrig,
        in_y1_2deg_irrig,
        in_y0_2deg_agri,
        in_y1_2deg_agri,
        out_y1_2deg_agri,
        np.ones_like(max_nfert_y1),
        debug_ij,
        "irrig",
        db_crop_idx,
    )
    check_bad_vals(None, mid_nfert, mid_irrig, None, [], "mid_y1_2deg", out_prec)

    if debug_nfert and debug_ij is not None:
        debug_out_deltas(
            "outy0_to_midy1",
            "Nfert",
            out_y0_2deg_nfert * out_y0_2deg_agri,
            mid_nfert * out_y1_2deg_agri,
            debug_ij,
            lpjg_crops,
        )
    if debug_nfert and debug_ij is not None and db_crop_idx is not None:
        i = debug_ij[0] - 1
        j = debug_ij[1] - 1
        print(
            "        DEBUG_NFERT_CELL:",
            f"ij=({debug_ij[0]},{debug_ij[1]})",
            f"crop={lpjg_crops[db_crop_idx]}",
            f"out_y0={out_y0_2deg_nfert[i,j,db_crop_idx]:.6e}",
            f"in_y0={in_y0_2deg_nfert[i,j,db_crop_idx]:.6e}",
            f"in_y1={in_y1_2deg_nfert[i,j,db_crop_idx]:.6e}",
            f"mid={mid_nfert[i,j,db_crop_idx]:.6e}",
            f"unm={unm_nfert[i,j,db_crop_idx]:.6e}",
            f"area_out_y1={out_y1_2deg_agri[i,j,db_crop_idx]:.6e}",
            f"area_out_y0={out_y0_2deg_agri[i,j,db_crop_idx]:.6e}",
        )
        max_rate = float(max_nfert_y1[db_crop_idx]) if max_nfert_y1 is not None else 0.0
        print(
            "        DEBUG_NFERT_CELL_FULL:",
            f"out_y0={float(out_y0_2deg_nfert[i,j,db_crop_idx]):.17e}",
            f"in_y0={float(in_y0_2deg_nfert[i,j,db_crop_idx]):.17e}",
            f"in_y1={float(in_y1_2deg_nfert[i,j,db_crop_idx]):.17e}",
            f"mid={float(mid_nfert[i,j,db_crop_idx]):.17e}",
            f"unm={float(unm_nfert[i,j,db_crop_idx]):.17e}",
            f"area_out_y1={float(out_y1_2deg_agri[i,j,db_crop_idx]):.17e}",
            f"area_out_y0={float(out_y0_2deg_agri[i,j,db_crop_idx]):.17e}",
            f"max_rate={max_rate:.17e}",
        )
    if debug_irrig and debug_ij is not None:
        debug_out_deltas(
            "outy0_to_midy1",
            "irrig",
            out_y0_2deg_irrig * out_y0_2deg_agri,
            mid_irrig * out_y1_2deg_agri,
            debug_ij,
            lpjg_crops,
        )
    if debug_irrig and debug_ij is not None and db_crop_idx is not None:
        i = debug_ij[0] - 1
        j = debug_ij[1] - 1
        print(
            "        DEBUG_IRRIG_CELL:",
            f"ij=({debug_ij[0]},{debug_ij[1]})",
            f"crop={lpjg_crops[db_crop_idx]}",
            f"out_y0={out_y0_2deg_irrig[i,j,db_crop_idx]:.6e}",
            f"in_y0={in_y0_2deg_irrig[i,j,db_crop_idx]:.6e}",
            f"in_y1={in_y1_2deg_irrig[i,j,db_crop_idx]:.6e}",
            f"mid={mid_irrig[i,j,db_crop_idx]:.6e}",
            f"unm={unm_irrig[i,j,db_crop_idx]:.6e}",
            f"area_out_y1={out_y1_2deg_agri[i,j,db_crop_idx]:.6e}",
            f"area_out_y0={out_y0_2deg_agri[i,j,db_crop_idx]:.6e}",
        )
        print(
            "        DEBUG_IRRIG_CELL_FULL:",
            f"out_y0={float(out_y0_2deg_irrig[i,j,db_crop_idx]):.17e}",
            f"in_y0={float(in_y0_2deg_irrig[i,j,db_crop_idx]):.17e}",
            f"in_y1={float(in_y1_2deg_irrig[i,j,db_crop_idx]):.17e}",
            f"mid={float(mid_irrig[i,j,db_crop_idx]):.17e}",
            f"unm={float(unm_irrig[i,j,db_crop_idx]):.17e}",
            f"area_out_y1={float(out_y1_2deg_agri[i,j,db_crop_idx]):.17e}",
            f"area_out_y0={float(out_y0_2deg_agri[i,j,db_crop_idx]):.17e}",
            "max_rate=1.00000000000000000e+00",
        )

    check_mgmt_conservation(
        out_y0_2deg_nfert,
        out_y0_2deg_agri,
        mid_nfert,
        out_y1_2deg_agri,
        in_y0_2deg_nfert,
        in_y0_2deg_agri,
        in_y1_2deg_nfert,
        in_y1_2deg_agri,
        unm_nfert,
        lpjg_crops,
        conserv_tol_pct,
        np.zeros(len(lpjg_crops), dtype=bool),
        "2b nfert",
        True,
    )
    check_mgmt_conservation(
        out_y0_2deg_irrig,
        out_y0_2deg_agri,
        mid_irrig,
        out_y1_2deg_agri,
        in_y0_2deg_irrig,
        in_y0_2deg_agri,
        in_y1_2deg_irrig,
        in_y1_2deg_agri,
        unm_irrig,
        lpjg_crops,
        conserv_tol_pct,
        np.zeros(len(lpjg_crops), dtype=bool),
        "2b irrig",
        True,
    )

    out_nfert, unm2_nfert, not_enough_nfert = ring_redistribute_mgmt(
        mid_nfert,
        out_y1_2deg_agri,
        unm_nfert,
        max_nfert_y1,
        lpjg_crops,
        debug_ij,
        {"maps_YXv": out_y0_2deg_nfert},
        out_y0_2deg_agri,
        {"maps_YXv": in_y0_2deg_nfert},
        in_y0_2deg_agri,
        {"maps_YXv": in_y1_2deg_nfert},
        in_y1_2deg_agri,
        conserv_tol_pct,
        "ringRedist nfert",
        db_crop_idx,
    )
    out_irrig, unm2_irrig, not_enough_irrig = ring_redistribute_mgmt(
        mid_irrig,
        out_y1_2deg_agri,
        unm_irrig,
        np.ones_like(max_nfert_y1),
        lpjg_crops,
        debug_ij,
        {"maps_YXv": out_y0_2deg_irrig},
        out_y0_2deg_agri,
        {"maps_YXv": in_y0_2deg_irrig},
        in_y0_2deg_agri,
        {"maps_YXv": in_y1_2deg_irrig},
        in_y1_2deg_agri,
        conserv_tol_pct,
        "ringRedist irrig",
        db_crop_idx,
    )
    check_bad_vals(None, out_nfert, out_irrig, None, [], "out_y1_2deg", out_prec)

    if debug_nfert and debug_ij is not None and db_crop_idx is not None:
        i = debug_ij[0] - 1
        j = debug_ij[1] - 1
        print(
            "        DEBUG_NFERT_CELL_OUT:",
            f"ij=({debug_ij[0]},{debug_ij[1]})",
            f"crop={lpjg_crops[db_crop_idx]}",
            f"out={out_nfert[i,j,db_crop_idx]:.6e}",
            f"unm2={unm2_nfert[i,j,db_crop_idx]:.6e}",
        )
        print(
            "        DEBUG_NFERT_CELL_OUT_FULL:",
            f"out={float(out_nfert[i,j,db_crop_idx]):.17e}",
            f"unm2={float(unm2_nfert[i,j,db_crop_idx]):.17e}",
        )
    if debug_irrig and debug_ij is not None and db_crop_idx is not None:
        i = debug_ij[0] - 1
        j = debug_ij[1] - 1
        print(
            "        DEBUG_IRRIG_CELL_OUT:",
            f"ij=({debug_ij[0]},{debug_ij[1]})",
            f"crop={lpjg_crops[db_crop_idx]}",
            f"out={out_irrig[i,j,db_crop_idx]:.6e}",
            f"unm2={unm2_irrig[i,j,db_crop_idx]:.6e}",
        )
        print(
            "        DEBUG_IRRIG_CELL_OUT_FULL:",
            f"out={float(out_irrig[i,j,db_crop_idx]):.17e}",
            f"unm2={float(unm2_irrig[i,j,db_crop_idx]):.17e}",
        )

    return (
        out_nfert,
        out_irrig,
        unm2_nfert,
        unm2_irrig,
        not_enough_nfert,
        not_enough_irrig,
        mid_nfert,
        mid_irrig,
        unm_nfert,
        unm_irrig,
    )


def get_unmet_mgmt(
    out_y0_2deg_mgmt_yxv: np.ndarray,
    out_y0_2deg_crop_yxv: np.ndarray,
    in_y0_2deg_mgmt_yxv: np.ndarray,
    in_y1_2deg_mgmt_yxv: np.ndarray,
    in_y0_2deg_crop_yxv: np.ndarray,
    in_y1_2deg_crop_yxv: np.ndarray,
    mid_y1_2deg_crop_yxv: np.ndarray,
    max_mgmt: np.ndarray,
    debug_ij: Tuple[int, int] | None = None,
    mgmt_name: str | None = None,
    db_crop: int | None = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute unmet mgmt (PLUMharm_getUnmet_mgmt.m)."""
    mgmt_tot_tol = 1e-6

    do_debug = debug_ij is not None and db_crop is not None
    if do_debug:
        dbi = debug_ij[0] - 1
        dbj = debug_ij[1] - 1
        print(
            f"PLUMharm_getUnmet_mgmt {mgmt_name}, [{debug_ij[0]} {debug_ij[1]}] LU {db_crop} (mgmtTot values):"
        )

    out_y0_tot = out_y0_2deg_mgmt_yxv * out_y0_2deg_crop_yxv
    in_y0_tot = in_y0_2deg_mgmt_yxv * in_y0_2deg_crop_yxv
    in_y1_tot = in_y1_2deg_mgmt_yxv * in_y1_2deg_crop_yxv
    
    if do_debug:
        print(f"   out_y0_tot:      {out_y0_tot[dbi, dbj, db_crop]:0.4e}")
        print(f"   in_y0_tot:       {in_y0_tot[dbi, dbj, db_crop]:0.4e}")
        print(f"   in_y1_tot:       {in_y1_tot[dbi, dbj, db_crop]:0.4e}")
    mgmt_d = in_y1_tot - in_y0_tot
    if do_debug:
        print(f"   Initial delta:    {mgmt_d[dbi, dbj, db_crop]:0.4e}")

    is_too_much_loss = out_y0_tot + mgmt_d < -mgmt_tot_tol
    any_too_much_loss = is_too_much_loss.copy()
    unmet_a = np.zeros_like(out_y0_tot)
    i = 0
    while np.any(is_too_much_loss):
        i += 1
        if i == 50:
            raise RuntimeError("Infinite loop?")
        too_much_by = (out_y0_tot + mgmt_d) * is_too_much_loss
        unmet_a = unmet_a + too_much_by
        mgmt_d = mgmt_d - too_much_by
        is_too_much_loss = out_y0_tot + mgmt_d < -mgmt_tot_tol
        any_too_much_loss |= is_too_much_loss

    mid1_tot = out_y0_tot + mgmt_d
    mid1_tot[mid1_tot < 0] = 0
    
    if do_debug:
        print(f"   After avoiding <0:   {mid1_tot[dbi, dbj, db_crop]:0.4e}")
        print(f"              UnmetA:   {unmet_a[dbi, dbj, db_crop]:0.4e}")

    max_mgmt_yxv = np.broadcast_to(max_mgmt, mid_y1_2deg_crop_yxv.shape)
    max_tot = max_mgmt_yxv * mid_y1_2deg_crop_yxv
    
    if do_debug:
        print(f"   Crop area:   {mid_y1_2deg_crop_yxv[dbi, dbj, db_crop]:0.4e}")
        print(f"   Max mgmt:    {max_tot[dbi, dbj, db_crop]:0.4e}")
        print(
            f"   Max mgmt rate: {max_mgmt_yxv[dbi, dbj, db_crop]:0.17e} "
            f"area_y1: {mid_y1_2deg_crop_yxv[dbi, dbj, db_crop]:0.17e}"
        )
    is_too_much_now = mid1_tot > max_tot
    unmet_b = (mid1_tot - max_tot) * is_too_much_now
    mid2_tot = mid1_tot.copy()
    mid2_tot[is_too_much_now] = max_tot[is_too_much_now]
    if np.any(mid2_tot > max_tot):
        raise RuntimeError("Go back to WHILE loop for upper-limiting mgmt inputs?")

    mid2_y1 = mid2_tot / mid_y1_2deg_crop_yxv
    mid2_y1[mid_y1_2deg_crop_yxv == 0] = 0
    
    mid2_y1[mid2_y1 > max_mgmt_yxv] = max_mgmt_yxv[mid2_y1 > max_mgmt_yxv]

    unmet = unmet_a + unmet_b
    if do_debug:
        print(f"   After avoiding >max: {mid2_tot[dbi, dbj, db_crop]:0.4e}")
        print(f"                UnmetB: {unmet_b[dbi, dbj, db_crop]:0.4e}")
        print(f"                 Unmet: {unmet[dbi, dbj, db_crop]:0.4e}")
        print(f"   mid2_y1:        {mid2_y1[dbi, dbj, db_crop]:0.4e}")

    if np.any(np.isnan(mid2_y1)):
        raise RuntimeError("Some NaN in mid2_y1_2deg_mgmt_YXv!")
    if np.any(np.isnan(unmet)):
        raise RuntimeError("Some NaN in unmetMgmt_YXv!")

    return unmet, mid2_y1
