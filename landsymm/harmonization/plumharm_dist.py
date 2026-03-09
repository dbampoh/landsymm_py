"""Distribution of deltas to finer resolution (MATLAB parity)."""
from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np


def distribute_area_deltas(
    land_area_yx: np.ndarray,
    land_area_2deg_yx: np.ndarray,
    out_y0_2deg_agri_yxv: np.ndarray,
    out_y1_2deg_agri_yxv: np.ndarray,
    out_y0_agri_yxv: np.ndarray,
    out_y0_vegd_yx: np.ndarray,
    fix_tiny_negs_tol_m2: float,
    conserv_tol_pct: float,
    conserv_tol_area: float,
    lu_names_agri: Sequence[str],
    debug_ijk: Tuple[int, int, int] | None = None,
) -> np.ndarray:
    """Distribute 2-deg deltas to 0.5-deg (PLUMharm_distDeltas_areaCrops_recursive.m)."""
    update_avail_land = True
    proper_zero_denoms = True

    do_debug = debug_ijk is not None
    if not do_debug:
        debug_ijk = (np.inf, np.inf, np.inf)

    n_agri = out_y0_agri_yxv.shape[2]
    out_y1_agri_yxv = out_y0_agri_yxv.copy()

    ny2, nx2 = land_area_2deg_yx.shape
    for i in range(ny2):
        for j in range(nx2):
            iy = np.arange(4 * i + 3, 4 * i - 1, -1)
            ix = np.arange(4 * j, 4 * j + 4)

            agri_d = out_y1_2deg_agri_yxv[i, j, :] - out_y0_2deg_agri_yxv[i, j, :]
            if np.all(agri_d == 0):
                out_y1_agri_yxv[np.ix_(iy, ix)] = out_y0_agri_yxv[np.ix_(iy, ix)]
                continue

            out_y0_vegd_these = out_y0_vegd_yx[np.ix_(iy, ix)].ravel(order="F")

            if update_avail_land:
                now_agri_yx = np.sum(out_y0_agri_yxv[np.ix_(iy, ix)], axis=2)
                avail_land = out_y0_vegd_these - now_agri_yx.ravel(order="F")
                if np.sum(avail_land) == 0 and np.sum(agri_d) > 0:
                    if np.sum(np.abs(agri_d)) < conserv_tol_area:
                        out_y1_agri_yxv[np.ix_(iy, ix)] = out_y0_agri_yxv[np.ix_(iy, ix)]
                        print(
                            f"Warning: ({i+1},{j+1}): Positive sum(agri_d) but no avail_land; "
                            f"sum(abs(agri_d)) < conserv_tol_area {conserv_tol_area}, setting agri_d to 0."
                        )
                        continue
                    print(
                        f"Warning: ({i+1},{j+1}): Positive sum(agri_d) but no avail_land; "
                        f"some abs(agri_d) >= conserv_tol_area {conserv_tol_area}."
                    )
            else:
                avail_land = np.array([])

            already_done = np.zeros(n_agri, dtype=bool)
            out_y1_agri_yxv, _ = _loop_thru_agri(
                already_done,
                0,
                debug_ijk,
                update_avail_land,
                proper_zero_denoms,
                conserv_tol_area,
                i,
                j,
                iy,
                ix,
                out_y0_vegd_these,
                agri_d,
                out_y0_agri_yxv,
                out_y0_2deg_agri_yxv,
                out_y1_agri_yxv,
                avail_land,
                lu_names_agri,
                fix_tiny_negs_tol_m2,
            )

            agri_yx = np.sum(out_y1_agri_yxv[np.ix_(iy, ix)], axis=2)
            max_too_much = np.max(agri_yx.ravel(order="F") - out_y0_vegd_these)
            if np.any(agri_yx.ravel(order="F") > fix_tiny_negs_tol_m2 + out_y0_vegd_these):
                print(
                    f"Warning: ({i+1},{j+1}): Members > vegd_area (by max {max_too_much}) "
                    "in half-deg out_y1_agri"
                )

            for k in range(n_agri):
                this_d_half = np.sum(
                    out_y1_agri_yxv[np.ix_(iy, ix, [k])] - out_y0_agri_yxv[np.ix_(iy, ix, [k])]
                )
                if agri_d[k] != 0:
                    if abs((this_d_half - agri_d[k]) / agri_d[k] * 100) > conserv_tol_pct and agri_d[
                        k
                    ] > 1e6:
                        raise RuntimeError(
                            f"Global {lu_names_agri[k]} area changes not conserved within "
                            f"{conserv_tol_pct}% (step 4)"
                        )

    return out_y1_agri_yxv


def distribute_mgmt(*args, **kwargs):
    """Distribute management to finer resolution (PLUMharm_distMgmt.m)."""
    in_yxv = args[0]
    in_res = args[1]
    out_res = args[2]
    res_ratio = in_res / out_res
    if not _is_power_of_two(res_ratio):
        raise RuntimeError("in_res / out_res must be a power of two!")
    res_ratio = int(res_ratio)
    out = in_yxv
    for _ in range(int(np.log2(res_ratio))):
        out = np.repeat(out, 2, axis=0)
        out = np.repeat(out, 2, axis=1)
    return out


def _is_power_of_two(val: float) -> bool:
    if val <= 0:
        return False
    log2 = np.log2(val)
    return float(int(log2)) == float(log2)


def _loop_thru_agri(
    already_done: np.ndarray,
    k1: int,
    debug_ijk: Tuple[int, int, int],
    update_avail_land: bool,
    proper_zero_denoms: bool,
    conserv_tol_area: float,
    i: int,
    j: int,
    iy: np.ndarray,
    ix: np.ndarray,
    out_y0_vegd_these: np.ndarray,
    agri_d: np.ndarray,
    out_y0_agri_yxv: np.ndarray,
    out_y0_2deg_agri_yxv: np.ndarray,
    out_y1_agri_yxv: np.ndarray,
    avail_land: np.ndarray,
    agri_names: Sequence[str],
    fix_tiny_negs_tol_m2: float,
) -> Tuple[np.ndarray, np.ndarray]:
    n_agri = out_y0_agri_yxv.shape[2]
    if k1 >= n_agri:
        raise RuntimeError("k1 > Nagri")

    for k in range(k1, n_agri):
        if update_avail_land:
            now_agri_yx = np.sum(out_y1_agri_yxv[np.ix_(iy, ix)], axis=2)
            avail_land = out_y0_vegd_these - now_agri_yx.ravel(order="F")

        if already_done[k]:
            out_y1_this = out_y1_agri_yxv[np.ix_(iy, ix, [k])].squeeze(axis=2)
        else:
            this_d = agri_d[k]
            if this_d == 0:
                out_y1_this = out_y0_agri_yxv[np.ix_(iy, ix, [k])].squeeze(axis=2)
                out_y1_agri_yxv[np.ix_(iy, ix, [k])] = out_y1_this[:, :, None]
                already_done[k] = True
                continue

            out_y0_this = out_y0_agri_yxv[np.ix_(iy, ix, [k])].squeeze(axis=2)
            out_y0_2deg_this = out_y0_2deg_agri_yxv[i, j, k]

            if this_d < 0:
                if proper_zero_denoms:
                    this_p = this_d / out_y0_2deg_this
                    if out_y0_2deg_this == 0:
                        this_p = 0
                else:
                    this_p = this_d / (out_y0_2deg_this + 1e-12)
                tmp = out_y0_this * (1 + this_p)
                out_y1_this = tmp
            elif this_d > 0:
                if np.sum(avail_land) < this_d and k < n_agri - 1:
                    if not update_avail_land:
                        now_agri_yx0 = np.sum(out_y0_agri_yxv[np.ix_(iy, ix)], axis=2)
                        avail_land_tmp = out_y0_vegd_these - now_agri_yx0.ravel(order="F")
                    else:
                        avail_land_tmp = np.array([])
                    out_y1_agri_yxv, already_done = _loop_thru_agri(
                        already_done,
                        k + 1,
                        debug_ijk,
                        update_avail_land,
                        proper_zero_denoms,
                        conserv_tol_area,
                        i,
                        j,
                        iy,
                        ix,
                        out_y0_vegd_these,
                        agri_d,
                        out_y0_agri_yxv,
                        out_y0_2deg_agri_yxv,
                        out_y1_agri_yxv,
                        avail_land_tmp,
                        agri_names,
                        fix_tiny_negs_tol_m2,
                    )
                    now_agri_yx_this = np.sum(out_y1_agri_yxv, axis=2)
                    avail_land_this = out_y0_vegd_these - now_agri_yx_this[np.ix_(iy, ix)].ravel(
                        order="F"
                    )
                    if np.sum(avail_land_this) + abs(fix_tiny_negs_tol_m2) < this_d:
                        print(
                            f"Warning: ({i+1},{j+1}): {agri_names[k]} increase {this_d} "
                            f"exceeds available land {np.sum(avail_land_this)} by "
                            f"{this_d - np.sum(avail_land_this)}"
                        )
                    if proper_zero_denoms:
                        if np.sum(avail_land_this) == 0:
                            raise RuntimeError(
                                f"No available land in cell ({i+1},{j+1},{k+1})"
                            )
                        tmp = out_y0_this.ravel(order="F") + avail_land_this / np.sum(
                            avail_land_this
                        ) * this_d
                    else:
                        tmp = out_y0_this.ravel(order="F") + avail_land_this / (
                            np.sum(avail_land_this) + 1e-12
                        ) * this_d
                    out_y1_this = tmp.reshape(out_y0_this.shape, order="F")
                else:
                    if proper_zero_denoms:
                        if np.sum(avail_land) == 0:
                            raise RuntimeError(
                                f"No available land in cell ({i+1},{j+1},{k+1})"
                            )
                        tmp = out_y0_this.ravel(order="F") + avail_land / np.sum(avail_land) * this_d
                    else:
                        tmp = out_y0_this.ravel(order="F") + avail_land / (
                            np.sum(avail_land) + 1e-12
                        ) * this_d
                    out_y1_this = tmp.reshape(out_y0_this.shape, order="F")
            else:
                raise RuntimeError("How is this_d not < or > 0??")

        tmp_vec = out_y1_this.ravel(order="F")
        if np.any(tmp_vec < -conserv_tol_area):
            raise RuntimeError("Negative members of half-deg out_y1_this_YXv!")
        if np.any(tmp_vec > conserv_tol_area + out_y0_vegd_these):
            raise RuntimeError("Members > vegd_area in half-deg out_y1_this_YX!")

        out_y1_agri_yxv[np.ix_(iy, ix, [k])] = out_y1_this[:, :, None]
        already_done[k] = True

    return out_y1_agri_yxv, already_done
