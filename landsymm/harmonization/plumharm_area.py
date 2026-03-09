"""Area harmonization logic (partial, MATLAB parity)."""
from __future__ import annotations

from typing import Tuple

import numpy as np

from .plumharm_checks import check_area_conservation, check_bad_vals
from .plumharm_debug import debug_out_deltas
from .plumharm_dist import distribute_area_deltas
from .plumharm_ring_redist import ring_redistribute_area


def harmonize_area(*args, **kwargs):
    """Harmonize land-use area changes (PLUMharm.m area section)."""
    (
        in_y0_2deg_agri_yxv,
        in_y1_2deg_agri_yxv,
        out_y0_2deg_agri_yxv,
        out_y0_2deg_vegd_yx,
        out_y0_2deg_ntrl_yx,
        base_2deg_bare_yx,
        land_area_2deg_yx,
        res_area_2deg_yx,
        lu_names,
        lu_names_agri,
        is_agri,
        set_area_to_zero,
        fix_tiny_negs_tol_m2,
        conserv_tol_pct,
        conserv_tol_area,
        debug_ij,
        db_crop,
    ) = args

    db_crop_idx = None
    if db_crop:
        try:
            db_crop_idx = int(db_crop)
        except (TypeError, ValueError):
            db_crop_idx = (
                lu_names_agri.index(db_crop) if db_crop in lu_names_agri else None
            )

    def _compose_full(agri_yxv, ntrl_yx, bare_yx):
        out = np.zeros((*agri_yxv.shape[:2], len(lu_names)))
        out[:, :, is_agri] = agri_yxv
        out[:, :, lu_names.index("NATURAL")] = ntrl_yx
        out[:, :, lu_names.index("BARREN")] = bare_yx
        return out

    agri_d = in_y1_2deg_agri_yxv - in_y0_2deg_agri_yxv
    for i, do_zero in enumerate(set_area_to_zero):
        if do_zero:
            agri_d[:, :, i] = -out_y0_2deg_agri_yxv[:, :, i]

    mid1_y1_2deg_agri = out_y0_2deg_agri_yxv + agri_d
    if debug_ij is not None:
        # Experimental (trace-only): collapse ULP-scale positive residues from
        # cancellation at zero-crossings in area handoff arithmetic.
        cross_zero_tol = 8 * np.finfo(float).eps * np.maximum(
            np.maximum(np.abs(in_y0_2deg_agri_yxv), np.abs(in_y1_2deg_agri_yxv)),
            1.0,
        )
        tiny_cancel = (
            (mid1_y1_2deg_agri > 0)
            & (mid1_y1_2deg_agri <= cross_zero_tol)
            & (out_y0_2deg_agri_yxv > 0)
            & (agri_d < 0)
            & (np.abs(out_y0_2deg_agri_yxv + agri_d) <= cross_zero_tol)
        )
        if np.any(tiny_cancel):
            mid1_y1_2deg_agri[tiny_cancel] = 0.0

    mid1_y1_2deg_ntrl = out_y0_2deg_ntrl_yx - np.sum(agri_d, axis=2)
    if debug_ij is not None:
        out_y0_bare = land_area_2deg_yx - out_y0_2deg_vegd_yx
        mid1_bare = land_area_2deg_yx - (np.sum(mid1_y1_2deg_agri, axis=2) + mid1_y1_2deg_ntrl)
        debug_out_deltas(
            "outy0_to_mid1y1",
            "areas",
            _compose_full(out_y0_2deg_agri_yxv, out_y0_2deg_ntrl_yx, out_y0_bare),
            _compose_full(mid1_y1_2deg_agri, mid1_y1_2deg_ntrl, mid1_bare),
            debug_ij,
            lu_names,
        )
        if db_crop_idx is not None:
            i = debug_ij[0] - 1
            j = debug_ij[1] - 1
            out_y0_dbg = float(out_y0_2deg_agri_yxv[i, j, db_crop_idx])
            in_y0_dbg = float(in_y0_2deg_agri_yxv[i, j, db_crop_idx])
            in_y1_dbg = float(in_y1_2deg_agri_yxv[i, j, db_crop_idx])
            mid1_dbg = float(mid1_y1_2deg_agri[i, j, db_crop_idx])
            delta_in_dbg = in_y1_dbg - in_y0_dbg
            mid1_calc_dbg = out_y0_dbg + delta_in_dbg
            print(
                "        DEBUG_AREA_CELL_PRE:",
                f"ij=({debug_ij[0]},{debug_ij[1]})",
                f"crop={lu_names_agri[db_crop_idx]}",
                f"out_y0={out_y0_dbg:.6e}",
                f"in_y0={in_y0_dbg:.6e}",
                f"in_y1={in_y1_dbg:.6e}",
                f"mid1={mid1_dbg:.6e}",
            )
            print(
                "        DEBUG_AREA_CELL_PRE_FULL:",
                f"out_y0={out_y0_dbg:.17e}",
                f"in_y0={in_y0_dbg:.17e}",
                f"in_y1={in_y1_dbg:.17e}",
                f"delta_in={delta_in_dbg:.17e}",
                f"mid1_calc={mid1_calc_dbg:.17e}",
                f"mid1={mid1_dbg:.17e}",
                f"mid1_minus_calc={(mid1_dbg - mid1_calc_dbg):.17e}",
            )

    unmet_agri, mid_y1_2deg_agri, mid_y1_2deg_ntrl = get_unmet_crop_area_res(
        mid1_y1_2deg_agri,
        out_y0_2deg_vegd_yx,
        res_area_2deg_yx,
        np.sum(out_y0_2deg_agri_yxv, axis=2),
        debug_ij,
        db_crop_idx,
    )
    mid_y1_2deg_vegd = np.sum(mid_y1_2deg_agri, axis=2) + mid_y1_2deg_ntrl
    mid_y1_2deg_bare = land_area_2deg_yx - mid_y1_2deg_vegd
    if debug_ij is not None:
        out_y0_bare = land_area_2deg_yx - out_y0_2deg_vegd_yx
        debug_out_deltas(
            "outy0_to_midy1.pre.fixTinyNegs",
            "areas",
            _compose_full(out_y0_2deg_agri_yxv, out_y0_2deg_ntrl_yx, out_y0_bare),
            _compose_full(mid_y1_2deg_agri, mid_y1_2deg_ntrl, mid_y1_2deg_bare),
            debug_ij,
            lu_names,
        )
        if db_crop_idx is not None:
            i = debug_ij[0] - 1
            j = debug_ij[1] - 1
            print(
                "        DEBUG_AREA_CELL_POST_UNMET:",
                f"ij=({debug_ij[0]},{debug_ij[1]})",
                f"crop={lu_names_agri[db_crop_idx]}",
                f"mid={mid_y1_2deg_agri[i,j,db_crop_idx]:.6e}",
            )

    tmp = np.concatenate(
        [mid_y1_2deg_agri, mid_y1_2deg_ntrl[:, :, None], mid_y1_2deg_bare[:, :, None]],
        axis=2,
    )
    tmp = fix_tiny_negs(
        tmp,
        np.repeat(land_area_2deg_yx[:, :, None], tmp.shape[2], axis=2),
        lu_names,
        6,
        fix_tiny_negs_tol_m2,
        conserv_tol_area,
        debug_ij,
    )
    mid_y1_2deg_agri = tmp[:, :, : -2]
    mid_y1_2deg_ntrl = tmp[:, :, -2]
    mid_y1_2deg_bare = tmp[:, :, -1]
    if debug_ij is not None:
        out_y0_bare = land_area_2deg_yx - out_y0_2deg_vegd_yx
        debug_out_deltas(
            "outy0_to_midy1.post.fixTinyNegs",
            "areas",
            _compose_full(out_y0_2deg_agri_yxv, out_y0_2deg_ntrl_yx, out_y0_bare),
            _compose_full(mid_y1_2deg_agri, mid_y1_2deg_ntrl, mid_y1_2deg_bare),
            debug_ij,
            lu_names,
        )
        debug_out_deltas(
            "outy0_to_midy1",
            "areas",
            _compose_full(out_y0_2deg_agri_yxv, out_y0_2deg_ntrl_yx, out_y0_bare),
            _compose_full(mid_y1_2deg_agri, mid_y1_2deg_ntrl, mid_y1_2deg_bare),
            debug_ij,
            lu_names,
        )

    check_area_conservation(
        out_y0_2deg_agri_yxv,
        mid_y1_2deg_agri,
        in_y0_2deg_agri_yxv,
        in_y1_2deg_agri_yxv,
        unmet_agri,
        lu_names_agri,
        conserv_tol_pct,
        conserv_tol_area,
        "2",
        True,
    )

    out_y1_2deg_agri, out_y1_2deg_ntrl = ring_redistribute_area(
        mid_y1_2deg_agri,
        unmet_agri,
        debug_ij,
        conserv_tol_pct,
        "2b areas",
        in_y0_2deg_agri_yxv,
        in_y1_2deg_agri_yxv,
        out_y0_2deg_agri_yxv,
        mid_y1_2deg_ntrl,
        res_area_2deg_yx,
        lu_names_agri,
    )
    if debug_ij is not None and db_crop_idx is not None:
        i = debug_ij[0] - 1
        j = debug_ij[1] - 1
        print(
            "        DEBUG_AREA_CELL_POST_RING:",
            f"ij=({debug_ij[0]},{debug_ij[1]})",
            f"crop={lu_names_agri[db_crop_idx]}",
            f"out={out_y1_2deg_agri[i,j,db_crop_idx]:.6e}",
        )
    out_y1_2deg_vegd = np.sum(out_y1_2deg_agri, axis=2) + out_y1_2deg_ntrl
    out_y1_2deg_bare = land_area_2deg_yx - out_y1_2deg_vegd

    tmp = np.concatenate(
        [out_y1_2deg_agri, out_y1_2deg_ntrl[:, :, None], out_y1_2deg_bare[:, :, None]],
        axis=2,
    )
    tmp = fix_tiny_negs(
        tmp,
        np.repeat(land_area_2deg_yx[:, :, None], tmp.shape[2], axis=2),
        lu_names,
        6,
        fix_tiny_negs_tol_m2,
        conserv_tol_area,
        debug_ij,
    )
    out_y1_2deg_agri = tmp[:, :, : -2]
    out_y1_2deg_ntrl = tmp[:, :, -2]
    out_y1_2deg_bare = tmp[:, :, -1]
    if debug_ij is not None and db_crop_idx is not None:
        i = debug_ij[0] - 1
        j = debug_ij[1] - 1
        print(
            "        DEBUG_AREA_CELL_POST_FIX:",
            f"ij=({debug_ij[0]},{debug_ij[1]})",
            f"crop={lu_names_agri[db_crop_idx]}",
            f"out={out_y1_2deg_agri[i,j,db_crop_idx]:.6e}",
        )

    check_bad_vals(
        _compose_full(out_y1_2deg_agri, out_y1_2deg_ntrl, out_y1_2deg_bare),
        None,
        None,
        land_area_2deg_yx,
        lu_names,
        "out_y1_2deg",
        6,
    )

    check_area_conservation(
        out_y0_2deg_agri_yxv,
        out_y1_2deg_agri,
        in_y0_2deg_agri_yxv,
        in_y1_2deg_agri_yxv,
        np.zeros_like(out_y0_2deg_agri_yxv),
        lu_names_agri,
        conserv_tol_pct,
        conserv_tol_area,
        "3",
        True,
    )

    return out_y1_2deg_agri, out_y1_2deg_ntrl, out_y1_2deg_bare


def get_unmet_crop_area_res(
    mid1_y1_2deg_agri_yxv: np.ndarray,
    base_2deg_vegd_yx: np.ndarray,
    res_area_2deg_yx: np.ndarray,
    sum_agri_y0_yx: np.ndarray,
    debug_ij: Tuple[int, int] | None = None,
    db_crop: int | None = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute unmet agri area from reserved-area constraints (PLUMharm_getUnmet_cropAreaRes.m)."""
    proper_zero_denoms = True
    n_agri = mid1_y1_2deg_agri_yxv.shape[2]

    unres_area = base_2deg_vegd_yx - res_area_2deg_yx
    if np.any(unres_area > base_2deg_vegd_yx):
        raise RuntimeError("How is unreserved area > base vegetated area?")
    max_agri_y1 = np.maximum(unres_area, sum_agri_y0_yx)
    max_agri_y1_yxv = np.repeat(max_agri_y1[:, :, None], n_agri, axis=2)

    mid_y1_2deg_agri_yxv = mid1_y1_2deg_agri_yxv.copy()

    unmet_a = (mid_y1_2deg_agri_yxv - max_agri_y1_yxv) * (
        mid_y1_2deg_agri_yxv > max_agri_y1_yxv
    )
    mid_y1_2deg_agri_yxv = mid_y1_2deg_agri_yxv - unmet_a

    unmet_b = mid_y1_2deg_agri_yxv * (mid_y1_2deg_agri_yxv < 0)
    mid_y1_2deg_agri_yxv = mid_y1_2deg_agri_yxv - unmet_b

    sum_agri_y1 = np.sum(mid_y1_2deg_agri_yxv, axis=2)
    is_ok = sum_agri_y1 <= max_agri_y1
    unmet_d = np.zeros_like(mid_y1_2deg_agri_yxv)
    j = 0
    while np.any(~is_ok):
        j += 1
        if j > 50:
            raise RuntimeError(
                'Possible infinite loop in "Check that sum(agri_y1) does not exceed max agri. area."'
            )
        exceed_land = (sum_agri_y1 - max_agri_y1) * (sum_agri_y1 > max_agri_y1)
        exceed_land_yxv = np.repeat(exceed_land[:, :, None], n_agri, axis=2)
        if proper_zero_denoms:
            weights = mid_y1_2deg_agri_yxv / np.repeat(sum_agri_y1[:, :, None], n_agri, axis=2)
            unmet_dtmp = weights * exceed_land_yxv
            unmet_dtmp[np.repeat(sum_agri_y1[:, :, None], n_agri, axis=2) == 0] = 0
        else:
            unmet_dtmp = (
                mid_y1_2deg_agri_yxv
                / (np.repeat(sum_agri_y1[:, :, None], n_agri, axis=2) + 1e-12)
                * exceed_land_yxv
            )
        unmet_d += unmet_dtmp
        mid_y1_2deg_agri_yxv -= unmet_dtmp
        sum_agri_y1 = np.sum(mid_y1_2deg_agri_yxv, axis=2)
        is_ok = sum_agri_y1 <= max_agri_y1

    mid_y1_2deg_ntrl_yx = base_2deg_vegd_yx - np.sum(mid_y1_2deg_agri_yxv, axis=2)
    unmet_yxv = unmet_a + unmet_b + unmet_d

    if debug_ij is not None:
        dbi = debug_ij[0] - 1
        dbj = debug_ij[1] - 1
        print("   After reserved-area checks (debug cell):")
        print(f"      agri_y0:\t{sum_agri_y0_yx[dbi, dbj]:0.4e}")
        print(f"      agri_y1_mid:\t{sum_agri_y1[dbi, dbj]:0.4e}")
        if db_crop is not None:
            print(
                f"      unmetA_this:\t{unmet_a[dbi, dbj, db_crop]:0.4e} "
                f"unmetB_this:\t{unmet_b[dbi, dbj, db_crop]:0.4e} "
                f"unmetD_this:\t{unmet_d[dbi, dbj, db_crop]:0.4e} "
                f"total_unmet_this:\t{unmet_yxv[dbi, dbj, db_crop]:0.4e}"
            )
        print(f"      ntrl_y1_mid:\t{mid_y1_2deg_ntrl_yx[dbi, dbj]:0.4e}")

    return unmet_yxv, mid_y1_2deg_agri_yxv, mid_y1_2deg_ntrl_yx


def fix_tiny_negs(
    in_yxv: np.ndarray,
    land_area_yxv: np.ndarray,
    lu_names,
    out_prec: int,
    fix_tiny_negs_tol_m2: float,
    conserv_tol_area: float,
    debug_ij: Tuple[int, int] | None = None,
) -> np.ndarray:
    """Fix tiny negative areas (PLUMharm_fixTinyNegs.m)."""
    n_lu = in_yxv.shape[2]
    fix_tiny_negs_tol_m2 = abs(fix_tiny_negs_tol_m2)
    conserv_tol_area = abs(conserv_tol_area)
    if fix_tiny_negs_tol_m2 > conserv_tol_area:
        raise RuntimeError("fixTinyNegs_tol_m2 must be <= conserv_tol_area")

    if np.min(in_yxv) < -fix_tiny_negs_tol_m2:
        msg_lines = [
            f"Negative value(s) of area exceeding tolerance {fix_tiny_negs_tol_m2} m2:"
        ]
        non_ntrl_bad = False
        for v in range(n_lu):
            this_min = np.min(in_yxv[:, :, v])
            if this_min >= -fix_tiny_negs_tol_m2:
                continue
            if len(lu_names) != n_lu:
                raise RuntimeError("length(LUnames) ~= Nlu")
            this_lu = lu_names[v]
            if this_lu != "NATURAL":
                non_ntrl_bad = True
            msg_lines.append(f"    {this_lu} min\t{this_min}")
        msg = "\n".join(msg_lines)
        if np.min(in_yxv) < -conserv_tol_area:
            raise RuntimeError(msg)
        msg = f"{msg}\nStill within conserv_tol_area ({conserv_tol_area})"
        if non_ntrl_bad:
            raise RuntimeError(
                f"{msg}, but this is only expected for NATURAL. Figure out why this happened!"
            )
        print(f"Warning: {msg}, and it's only for NATURAL, so zeroing.")

    out = in_yxv.copy()
    is_bad = np.any(out < 0, axis=2)
    if np.any(is_bad):
        n_bad = np.sum(is_bad)
        is_bad_yxv = np.repeat(is_bad[:, :, None], n_lu, axis=2)
        tmp = out[is_bad_yxv]
        tmp_land = land_area_yxv[is_bad_yxv]
        tmp = tmp / tmp_land
        if np.any(np.isnan(tmp)):
            raise RuntimeError("How do you have NaN in tmp?")
        tmp_xv = tmp.reshape((n_bad, n_lu))
        tmp_xv[tmp_xv < 0] = 0
        tmp_sum = np.sum(tmp_xv, axis=1)
        j = 0
        while np.max(np.abs(tmp_sum - 1)) > 3 * np.finfo(float).eps:
            j += 1
            if j > 50:
                raise RuntimeError("Possible infinite loop in fixing tiny negative areas!")
            tmp_xv = tmp_xv / tmp_sum[:, None]
            tmp_sum = np.sum(tmp_xv, axis=1)
        out[is_bad_yxv] = tmp_xv.ravel() * tmp_land
        small_diff = np.abs(out - in_yxv) < 1e-6
        out[small_diff] = in_yxv[small_diff]

    small_neg = (out < 0) & (out > -out_prec / 2)
    if np.any(small_neg):
        total_negs = np.sum(out[small_neg])
        if abs(total_negs) >= 1:
            print(f"Warning: Zeroing out {total_negs:0.1g} total small negatives")
        out[small_neg] = 0
    return out
