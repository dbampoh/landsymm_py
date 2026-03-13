"""Conservation and validation checks (MATLAB parity)."""
from __future__ import annotations

from typing import Sequence

import numpy as np


def check_area_conservation(
    out_y0_area,
    out_y1_area,
    in_y0_area,
    in_y1_area,
    total_unmet_area,
    lu_names: Sequence[str],
    tol_pct,
    tol_area,
    check_name,
    do_warn,
):
    """Check area conservation (PLUMharm_checkCons_area.m)."""
    area_d = in_y1_area - in_y0_area
    bad = 0
    n_lu = len(lu_names)
    for i in range(n_lu):
        area_d2 = (out_y1_area[:, :, i] + total_unmet_area[:, :, i]) - out_y0_area[:, :, i]
        area_d_glob_1 = np.sum(area_d[:, :, i])
        area_d_glob_2 = np.sum(area_d2)
        if np.isnan(area_d_glob_1):
            raise RuntimeError(f"NaN produced in check {check_name}, area_d_glob_1")
        if np.isnan(area_d_glob_2):
            raise RuntimeError(f"NaN produced in check {check_name}, area_d_glob_2")
        err = area_d_glob_2 - area_d_glob_1
        if area_d_glob_1 == 0:
            local_tol_area = 1
            is_bad = abs(err) > local_tol_area
            err_pct = np.inf
        else:
            err_pct = err / area_d_glob_1 * 100
            is_bad = abs(err_pct) > tol_pct
        if is_bad and abs(err_pct) > tol_area:
            if np.sum(out_y1_area[:, :, i]) == 0:
                bad = -1
                if do_warn:
                    print(
                        "Warning: Global {} area changes are not conserved to within {} percent because they hit zero "
                        "(check {})".format(lu_names[i], f"{tol_pct:g}", check_name)
                    )
            else:
                bad = 1
                if do_warn:
                    if area_d_glob_1 == 0:
                        print(
                            "Warning: Global {} area changes are not conserved to within {} m2!\n"
                            "(check {}; err = {:0.2f})".format(
                                lu_names[i], f"{tol_area:g}", check_name, err
                            )
                        )
                    else:
                        print(
                            "Warning: Global {} area changes are not conserved to within {:0.2f} percent!\n"
                            "(check {}; err = {:0.2f})".format(
                                lu_names[i], tol_pct, check_name, err_pct
                            )
                        )
    return bad


def check_mgmt_conservation(
    out_y0_mgmt, out_y0_area, out_y1_mgmt, out_y1_area,
    in_y0_mgmt, in_y0_area, in_y1_mgmt, in_y1_area,
    total_unmet, crops, tol_pct, not_enough, check_name, do_warn
):
    """Check management conservation (PLUMharm_checkCons_mgmt.m)."""
    in_y0_tot = in_y0_mgmt * in_y0_area
    in_y1_tot = in_y1_mgmt * in_y1_area
    out_y0_tot = out_y0_mgmt * out_y0_area
    out_y1_tot = out_y1_mgmt * out_y1_area
    mgmt_d = in_y1_tot - in_y0_tot
    bad = 0
    for i, crop in enumerate(crops):
        mgmt_d2 = out_y1_tot[:, :, i] - out_y0_tot[:, :, i] + total_unmet[:, :, i]
        mgmt_d_glob_1 = np.sum(mgmt_d[:, :, i])
        mgmt_d_glob_2 = np.sum(mgmt_d2)
        if np.isnan(mgmt_d_glob_1):
            raise RuntimeError(f"NaN produced in check {check_name}, mgmt_d_glob_1")
        if np.isnan(mgmt_d_glob_2):
            raise RuntimeError(f"NaN produced in check {check_name}, mgmt_d_glob_2")
        if mgmt_d_glob_1 == 0:
            if mgmt_d_glob_2 != 0:
                if do_warn:
                    print(
                        "Warning: Global {} mgmt changes are not conserved because baseline delta is 0 "
                        "(check {}).".format(crop, check_name)
                    )
                bad = -1
            continue
        pct_error = (mgmt_d_glob_2 - mgmt_d_glob_1) / mgmt_d_glob_1 * 100
        if abs(pct_error) > tol_pct:
            if not_enough[i]:
                bad = -1
                if do_warn and np.sum(out_y1_tot[:, :, i]) == 0:
                    print(
                        "Warning: Global {} mgmt changes are not conserved to within {:0.2f} percent because "
                        "they hit zero (check {})".format(crop, tol_pct, check_name)
                    )
                elif do_warn and not_enough[i]:
                    print(
                        "Warning: Global {} mgmt changes are not conserved to within {:0.2f} percent because "
                        "not enough available headroom (check {})".format(crop, tol_pct, check_name)
                    )
            else:
                bad = 1
                raise RuntimeError(
                    "Global {} mgmt changes are not conserved to within {:0.2f}%! "
                    "(check {}; err = {:0.2f}%)".format(crop, tol_pct, check_name, pct_error)
                )
    return bad


def check_preserved_deltas(
    out_y0, out_y1, in_y0, in_y1, crops, check_name, do_warn, conserv_tol_pct=0.2
):
    """Check preserved global deltas with percentage tolerance (matching MATLAB)."""
    in_d = in_y1 - in_y0
    out_d = out_y1 - out_y0
    for i, crop in enumerate(crops):
        in_glob = np.sum(in_d[:, :, i])
        out_glob = np.sum(out_d[:, :, i])
        if in_glob == 0:
            continue
        diff_pct = 100 * (out_glob - in_glob) / abs(in_glob)
        if abs(diff_pct) > conserv_tol_pct and do_warn:
            print(
                f"Warning: Global {crop} delta not preserved (check {check_name}; "
                f"in={in_glob:0.6g}, out={out_glob:0.6g}, diff={diff_pct:0.2f}%)"
            )


def check_bad_vals(
    in_lu_yxv,
    in_nfert_yxv,
    in_irrig_yxv,
    land_area_yx,
    lu_names,
    msg,
    out_prec,
    allow_unveg: bool = False,
):
    """Mirror PLUMharm_checkBadVals.m."""
    if in_lu_yxv is not None:
        is_bare = np.array([n == "BARREN" for n in lu_names], dtype=bool)
        if np.any(is_bare):
            vegd_yx = np.sum(in_lu_yxv[:, :, ~is_bare], axis=2)
        if land_area_yx is None:
            raise RuntimeError("land_area_yx required for LU checks.")
        with np.errstate(divide="ignore", invalid="ignore"):
            in_lu_frac = in_lu_yxv / land_area_yx[:, :, None]
        in_lu_frac = np.where(land_area_yx[:, :, None] == 0, 0, in_lu_frac)
        if np.any(np.isnan(in_lu_yxv)):
            raise RuntimeError(f"NaN(s) in {msg} LU maps!")
        if np.min(in_lu_frac) < 0:
            raise RuntimeError(f"Negative value(s) in {msg} LU maps! Min {np.min(in_lu_frac):0.1e}")
        if np.max(in_lu_frac) > 1 + 10 ** (-out_prec):
            raise RuntimeError(
                f"Value(s) >1 in {msg} LU maps! Max overage {np.max(in_lu_frac) - 1:0.1e}"
            )
        if np.max(np.sum(in_lu_frac, axis=2)) > 1 + 10 ** (-out_prec):
            raise RuntimeError(
                f"Sum(s) >1 in {msg} LU maps! Max overage {np.max(np.sum(in_lu_frac, axis=2)) - 1:0.1e}"
            )
        if (
            not allow_unveg
            and np.any(is_bare)
            and np.any((vegd_yx == 0) & (land_area_yx > 0))
        ):
            max_land = np.max(land_area_yx[vegd_yx == 0])
            tot_land = np.sum(land_area_yx[vegd_yx == 0])
            raise RuntimeError(
                f"Zero vegetation in cell(s) with land! Max land area {max_land:0.1e}, total {tot_land:0.1e}."
            )

    if in_nfert_yxv is not None:
        if np.any(np.isnan(in_nfert_yxv)):
            raise RuntimeError(f"NaN(s) in {msg} nfert maps!")
        if np.any(in_nfert_yxv < 0):
            raise RuntimeError(
                f"Negative value(s) in {msg} nfert maps! Min {np.min(in_nfert_yxv):0.1e}"
            )

    if in_irrig_yxv is not None:
        if np.any(np.isnan(in_irrig_yxv)):
            raise RuntimeError(f"NaN(s) in {msg} irrig maps!")
        if np.any(in_irrig_yxv < 0):
            raise RuntimeError(
                f"Negative value(s) in {msg} irrig maps! Min {np.min(in_irrig_yxv):0.1e}"
            )
        if np.max(in_irrig_yxv) > 1 + 10**0:
            raise RuntimeError(
                f"Value(s) >1 in {msg} irrig maps! Max overage {np.max(in_irrig_yxv) - 1:0.1e}"
            )
