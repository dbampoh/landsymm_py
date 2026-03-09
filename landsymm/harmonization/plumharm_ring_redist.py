"""Ring redistribution helpers (MATLAB parity)."""
from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np


def _sum_f(arr: np.ndarray) -> float:
    """Sum in Fortran order to match MATLAB reduction order."""
    total = 0.0
    for val in arr.ravel(order="F"):
        total += float(val)
    return total


def ring_redistribute_area(
    mid_y1_2deg_agri_yxv: np.ndarray,
    total_unmet_agri_yxv: np.ndarray,
    debug_ij: Tuple[int, int] | None,
    conserv_tol_pct: float,
    check_name: str,
    in_y0_area_yxv: np.ndarray,
    in_y1_area_yxv: np.ndarray,
    out_y0_2deg_agri_yxv: np.ndarray,
    out_y0_2deg_ntrl_yx: np.ndarray,
    res_area_2deg_yx: np.ndarray,
    lu_names_agri: Sequence[str],
) -> Tuple[np.ndarray, np.ndarray]:
    """Redistribute unmet agri area via rings (PLUMharm_ringRedist_areaCropsRes.m)."""
    out_ntrl_yx = out_y0_2deg_ntrl_yx.copy()
    out_y1_2deg_agri_yxv = mid_y1_2deg_agri_yxv.copy()

    ny, nx, n_agri = out_y1_2deg_agri_yxv.shape
    do_debug = debug_ij is not None
    if do_debug:
        displaced_agri_yxv = np.zeros_like(out_y1_2deg_agri_yxv)
        mean_dist_yxv = np.zeros_like(out_y1_2deg_agri_yxv)
        this_cell_of_int = (debug_ij[0] - 1) * nx + (debug_ij[1] - 1)
    else:
        displaced_agri_yxv = None
        mean_dist_yxv = None
        this_cell_of_int = None

    is_done = np.zeros_like(out_y1_2deg_agri_yxv, dtype=bool)

    loops = 0

    while np.any(~is_done):
        loops += 1
        if loops > 100:
            raise RuntimeError(f"Possible infinite loop (top while, {check_name}).")
        is_done_begin = is_done.copy()

        for k in range(ny):
            for m in range(nx):
                this_cell = (k * nx) + m
                nonres_ntrl = max(out_ntrl_yx.ravel()[this_cell] - res_area_2deg_yx.ravel()[this_cell], 0)
                avail_this = nonres_ntrl - np.sum(out_y1_2deg_agri_yxv[k, m, :])
                unmet_this = total_unmet_agri_yxv[k, m, :]

                if avail_this > 0 and np.any(unmet_this > 0):
                    unmet_pos = unmet_this * (unmet_this > 0)
                    unmet_red_tot = min(avail_this, np.sum(unmet_pos))
                    if np.sum(unmet_pos) > 0:
                        unmet_red = unmet_red_tot * unmet_pos / np.sum(unmet_pos)
                    else:
                        unmet_red = unmet_pos
                    out_y1_2deg_agri_yxv[k, m, :] += unmet_red
                    total_unmet_agri_yxv[k, m, :] -= unmet_red
                    out_ntrl_yx.ravel()[this_cell] -= np.sum(unmet_red)
                    if out_ntrl_yx.ravel()[this_cell] < 0:
                        raise RuntimeError("out_ntrl_yx(thisCell)<0")

                for i in range(n_agri):
                    if is_done[k, m, i]:
                        continue
                    if total_unmet_agri_yxv[k, m, i] == 0:
                        is_done[k, m, i] = True
                        continue

                    j = 0
                    out_this_yx = out_y1_2deg_agri_yxv[:, :, i].copy()
                    total_unmet_cell = total_unmet_agri_yxv[k, m, i]

                    while abs(total_unmet_cell) > 1e-8:
                        j += 1
                        if j > 100:
                            raise RuntimeError(
                                f"Possible infinite loop in crop ring adjustments: [{k+1} {m+1} {i+1}]"
                            )

                        out_agri_yx = np.sum(out_y1_2deg_agri_yxv, axis=2)
                        i_k = (k + np.arange(-j, j + 1)) % ny
                        i_m = (m + np.arange(-j, j + 1)) % nx
                        I_K, I_M = np.meshgrid(i_k, i_m, indexing="ij")
                        ring_lin = I_K.ravel(order="F") + I_M.ravel(order="F") * ny
                        if j * 2 + 1 > min(ny, nx):
                            ring_lin = np.unique(ring_lin)
                        ring_rows = ring_lin % ny
                        ring_cols = ring_lin // ny

                        if total_unmet_cell < 0:
                            avail_ring = out_this_yx[ring_rows, ring_cols]
                            total_avail = _sum_f(avail_ring)
                        elif total_unmet_cell > 0:
                            nonres_ring = np.maximum(
                                out_ntrl_yx[ring_rows, ring_cols]
                                - res_area_2deg_yx[ring_rows, ring_cols],
                                0,
                            )
                            avail_ring = nonres_ring - out_agri_yx[ring_rows, ring_cols]
                            avail_ring[avail_ring < 0] = 0
                            total_avail = _sum_f(avail_ring)
                        else:
                            raise RuntimeError("How was this not skipped?")

                        if total_avail == 0:
                            continue

                        out_this_yx, out_ntrl_yx, total_unmet_cell, displaced, mean_dist = _do_rings_area(
                            out_this_yx,
                            out_agri_yx,
                            total_unmet_cell,
                            displaced_agri_yxv[:, :, i] if do_debug else None,
                            (ring_rows, ring_cols),
                            out_ntrl_yx,
                            res_area_2deg_yx,
                            mean_dist_yxv[:, :, i] if do_debug else None,
                            this_cell_of_int,
                        )
                        out_y1_2deg_agri_yxv[:, :, i] = out_this_yx
                        if do_debug:
                            displaced_agri_yxv[:, :, i] = displaced
                            mean_dist_yxv[:, :, i] = mean_dist
                        total_unmet_agri_yxv[k, m, i] = total_unmet_cell

                        if do_debug:
                            from .plumharm_checks import \
                                check_area_conservation

                            bad = check_area_conservation(
                                out_y0_2deg_agri_yxv,
                                out_y1_2deg_agri_yxv,
                                in_y0_area_yxv,
                                in_y1_area_yxv,
                                total_unmet_agri_yxv,
                                lu_names_agri,
                                conserv_tol_pct,
                                0,
                                check_name,
                                True,
                            )
                            if bad == 1:
                                raise RuntimeError("Conservation check failed during debug.")

                    is_done[k, m, i] = True

        if np.array_equal(is_done_begin, is_done):
            raise RuntimeError("Not enough mgmt in the world!")

    return out_y1_2deg_agri_yxv, out_ntrl_yx


def ring_redistribute_mgmt(*args, **kwargs):
    """Redistribute unmet management via rings (TODO)."""
    (
        mid_y1_2deg_mgmt_yxv,
        out_y1_2deg_crop_area_yxv,
        total_unmet_mgmt_yxv,
        max_mgmt,
        lpjg_crops,
        debug_ij,
        out_y0_mgmt,
        out_y0_area_yxv,
        in_y0_mgmt,
        in_y0_area_yxv,
        in_y1_mgmt,
        in_y1_area_yxv,
        conserv_tol_pct,
        check_name,
        db_crop,
    ) = args

    if mid_y1_2deg_mgmt_yxv.shape != out_y1_2deg_crop_area_yxv.shape:
        raise RuntimeError("mid_y1_2deg_mgmt_yxv and out_y1_2deg_crop_area_yxv must be same size!")

    ny, nx, n_crops = mid_y1_2deg_mgmt_yxv.shape
    out_y1_2deg_mgmt_yxv = mid_y1_2deg_mgmt_yxv.copy()
    is_done = np.zeros_like(out_y1_2deg_mgmt_yxv, dtype=bool)
    not_enough = np.zeros(n_crops, dtype=bool)

    do_debug = debug_ij is not None and db_crop is not None
    dbk = debug_ij[0] - 1 if do_debug else None
    dbm = debug_ij[1] - 1 if do_debug else None
    if do_debug:
        print(
            f"ring_redistribute_mgmt debug cell: ij=({debug_ij[0]},{debug_ij[1]}) crop={db_crop}"
        )
        print(
            f"        debug total_unmet at cell =\t{total_unmet_mgmt_yxv[dbk, dbm, db_crop]:0.4e}"
        )
        tracked_val = out_y1_2deg_mgmt_yxv[dbk, dbm, db_crop]
    else:
        tracked_val = None

    loops = 0
    while np.any(~is_done):
        loops += 1
        if loops > 100:
            raise RuntimeError(f"Possible infinite loop (top while, {check_name}).")
        is_done_begin = is_done.copy()
        any_skipped = np.zeros(n_crops, dtype=bool)
        check_too_much = np.ones(n_crops, dtype=bool)

        for k in range(ny):
            for m in range(nx):
                this_cell = (k * nx) + m

                for i in range(n_crops):
                    if is_done[k, m, i]:
                        continue
                    if do_debug and k == dbk and m == dbm and i == db_crop:
                        print(
                            f"{lpjg_crops[i]}, pre-skip unmet =\t\t{total_unmet_mgmt_yxv[k, m, i]:0.4e}"
                        )

                    if check_too_much[i]:
                        check_too_much[i] = False
                        out_tot = out_y1_2deg_mgmt_yxv[:, :, i] * out_y1_2deg_crop_area_yxv[:, :, i]
                        max_tot = max_mgmt[i] * out_y1_2deg_crop_area_yxv[:, :, i]
                        total_unmet = total_unmet_mgmt_yxv[:, :, i]
                        total_pos_unmet = _sum_f(total_unmet[total_unmet > 0])
                        total_neg_unmet = _sum_f(total_unmet[total_unmet < 0])
                        total_avail_take = _sum_f(max_tot - out_tot) - total_neg_unmet
                        total_avail_give = _sum_f(out_tot) + total_pos_unmet
                        if do_debug and i == db_crop:
                            print(
                                f"check_too_much {lpjg_crops[i]}: pos_unmet={total_pos_unmet:0.4e} "
                                f"neg_unmet={total_neg_unmet:0.4e} avail_take={total_avail_take:0.4e} "
                                f"avail_give={total_avail_give:0.4e}"
                            )
                        too_much_pos = total_pos_unmet > total_avail_take
                        too_much_neg = -total_neg_unmet > total_avail_give
                        if too_much_pos and too_much_neg:
                            raise RuntimeError(
                                "Too much positive AND negative unmet; unhandled case."
                            )
                        if too_much_pos:
                            out_y1_2deg_mgmt_yxv[:, :, i] = np.where(
                                out_y1_2deg_crop_area_yxv[:, :, i] == 0,
                                0,
                                max_tot / out_y1_2deg_crop_area_yxv[:, :, i],
                            )
                            total_unmet_mgmt_yxv[:, :, i] = 0
                            is_done[:, :, i] = True
                            not_enough[i] = True
                            continue
                        if too_much_neg:
                            out_y1_2deg_mgmt_yxv[:, :, i] = 0
                            total_unmet_mgmt_yxv[:, :, i] = 0
                            is_done[:, :, i] = True
                            not_enough[i] = True
                            continue

                    max_tot_this = max_mgmt[i] * out_y1_2deg_crop_area_yxv[k, m, i]
                    out_total_this = out_y1_2deg_mgmt_yxv[k, m, i] * out_y1_2deg_crop_area_yxv[k, m, i]
                    if total_unmet_mgmt_yxv[k, m, i] > 0 and out_total_this < max_tot_this:
                        unmet_red = min(
                            max_tot_this - out_total_this, total_unmet_mgmt_yxv[k, m, i]
                        )
                        if out_y1_2deg_crop_area_yxv[k, m, i] > 0:
                            out_y1_2deg_mgmt_yxv[k, m, i] = (
                                out_total_this + unmet_red
                            ) / out_y1_2deg_crop_area_yxv[k, m, i]
                        total_unmet_mgmt_yxv[k, m, i] -= unmet_red

                    if total_unmet_mgmt_yxv[k, m, i] == 0:
                        if do_debug and k == dbk and m == dbm and i == db_crop:
                            print("                unmet is zero; skipping")
                        is_done[k, m, i] = True
                        continue

                    out_tot = out_y1_2deg_mgmt_yxv[:, :, i] * out_y1_2deg_crop_area_yxv[:, :, i]
                    if total_unmet_mgmt_yxv[k, m, i] < 0:
                        total_avail_world = _sum_f(out_tot)
                        if do_debug and k == dbk and m == dbm and i == db_crop:
                            print(
                                f"                total_avail_world =\t\t{total_avail_world:0.4e}"
                            )
                        if -total_unmet_mgmt_yxv[k, m, i] > total_avail_world:
                            if do_debug and k == dbk and m == dbm and i == db_crop:
                                print("                skipping: not enough mgmt in world (take)")
                            any_skipped[i] = True
                            continue
                    elif total_unmet_mgmt_yxv[k, m, i] > 0:
                        max_tot = max_mgmt[i] * out_y1_2deg_crop_area_yxv[:, :, i]
                        total_avail_world = _sum_f(max_tot - out_tot)
                        if do_debug and k == dbk and m == dbm and i == db_crop:
                            print(
                                f"                total_avail_world =\t\t{total_avail_world:0.4e}"
                            )
                        if total_unmet_mgmt_yxv[k, m, i] > total_avail_world:
                            if do_debug and k == dbk and m == dbm and i == db_crop:
                                print("                skipping: not enough headroom in world (give)")
                            any_skipped[i] = True
                            continue
                    else:
                        raise RuntimeError("If total_unmet == 0, you should have skipped!")

                    total_unmet_cell = total_unmet_mgmt_yxv[k, m, i]
                    if do_debug and k == dbk and m == dbm and i == db_crop:
                        print(
                            f"{lpjg_crops[i]}, j = 0, total_unmet_mgmt_YXv(k+1,m+1,i) =\t{total_unmet_cell:0.4e}"
                        )
                    j = 0
                    out_this_yx = out_y1_2deg_mgmt_yxv[:, :, i].copy()
                    out_this_area = out_y1_2deg_crop_area_yxv[:, :, i]
                    while abs(total_unmet_cell) > 1e-8:
                        j += 1
                        if j > 100:
                            raise RuntimeError(f"Possible infinite loop in mgmt ring adjustments ({check_name}).")

                        i_k = (k + np.arange(-j, j + 1)) % ny
                        i_m = (m + np.arange(-j, j + 1)) % nx
                        I_K, I_M = np.meshgrid(i_k, i_m, indexing="ij")
                        ring_lin = I_K.ravel(order="F") + I_M.ravel(order="F") * ny
                        if j * 2 + 1 > min(ny, nx):
                            ring_lin = np.unique(ring_lin)
                        ring_rows = ring_lin % ny
                        ring_cols = ring_lin // ny

                        out_total_ring = out_this_yx[ring_rows, ring_cols] * out_this_area[
                            ring_rows, ring_cols
                        ]
                        if total_unmet_cell < 0:
                            avail_ring = out_total_ring
                            total_avail_ring = _sum_f(avail_ring)
                        else:
                            max_tot_ring = max_mgmt[i] * out_this_area[ring_rows, ring_cols]
                            avail_ring = max_tot_ring - out_total_ring
                            avail_ring[avail_ring < 0] = 0
                            total_avail_ring = _sum_f(avail_ring)
                        if do_debug and k == dbk and m == dbm and i == db_crop:
                            if total_unmet_cell < 0:
                                total_avail_world = _sum_f(
                                    out_this_yx * out_this_area
                                )
                            else:
                                total_avail_world = _sum_f(
                                    max_mgmt[i] * out_this_area - out_this_yx * out_this_area
                                )
                            print(f"                                j =\t\t{j}")
                            print(
                                f"                total_avail_world =\t\t{total_avail_world:0.4e}"
                            )
                            print(
                                f"{lpjg_crops[i]}, j = {j}, total_unmet_mgmt_YXv(k+1,m+1,i) =\t{total_unmet_cell:0.4e}"
                            )
                            print(f"                total_avail_ring =\t\t{total_avail_ring:0.4e}")
                        if total_avail_ring == 0:
                            continue

                        if do_debug and k == dbk and m == dbm and i == db_crop:
                            cell_total = out_this_yx[k, m] * out_this_area[k, m]
                            print(
                                f"                pre do_rings: out_cell={out_this_yx[k, m]:0.4e} "
                                f"area_cell={out_this_area[k, m]:0.4e} total_cell={cell_total:0.4e}"
                            )

                        out_this_yx, total_unmet_cell = _do_rings_mgmt(
                            out_this_yx,
                            total_unmet_cell,
                            out_this_area,
                            max_mgmt[i],
                            (ring_rows, ring_cols),
                            debug=(do_debug and k == dbk and m == dbm and i == db_crop),
                            debug_prefix="                ",
                        )
                        if do_debug and k == dbk and m == dbm and i == db_crop:
                            print(
                                f"                after do_rings: out_cell={out_this_yx[k, m]:0.4e} "
                                f"unmet={total_unmet_cell:0.4e}"
                            )
                        out_y1_2deg_mgmt_yxv[:, :, i] = out_this_yx

                    is_done[k, m, i] = True
                    total_unmet_mgmt_yxv[k, m, i] = total_unmet_cell
                if do_debug:
                    current_val = out_y1_2deg_mgmt_yxv[dbk, dbm, db_crop]
                    if current_val != tracked_val:
                        print(
                            f"        debug cell updated by k={k+1}, m={m+1}: "
                            f"{tracked_val:0.6e} -> {current_val:0.6e}"
                        )
                        tracked_val = current_val

        if np.array_equal(is_done_begin, is_done):
            raise RuntimeError("Not enough mgmt in the world!")

    return out_y1_2deg_mgmt_yxv, total_unmet_mgmt_yxv, not_enough


def _do_rings_area(
    in_yx: np.ndarray,
    out_agri_yx: np.ndarray,
    total_unmet_in: float,
    displaced_yx: np.ndarray | None,
    ring_rc: Tuple[np.ndarray, np.ndarray],
    in_ntrl_yx: np.ndarray,
    res_area_yx: np.ndarray,
    mean_dist_yx: np.ndarray | None,
    this_cell_of_int: int | None,
):
    nonres_ntrl = np.maximum(in_ntrl_yx - res_area_yx, 0)
    out_yx = in_yx.copy()
    total_unmet_out = total_unmet_in

    ring_rows, ring_cols = ring_rc

    if total_unmet_in > 0:
        avail_space = nonres_ntrl[ring_rows, ring_cols] - out_agri_yx[ring_rows, ring_cols]
        avail_space[avail_space < 0] = 0
        total_avail = _sum_f(avail_space)
        if total_avail > 0:
            if total_unmet_in >= total_avail:
                to_ring = avail_space
                total_unmet_out = total_unmet_in - total_avail
            else:
                to_ring = total_unmet_in * avail_space / total_avail
                total_unmet_out = 0
            out_this_new = out_yx[ring_rows, ring_cols] + to_ring
            out_yx[ring_rows, ring_cols] = out_this_new
            if displaced_yx is not None:
                displaced_yx[ring_rows, ring_cols] += to_ring
                j_ring = (np.sqrt(len(ring_rows)) - 1) / 2
                now_weighted = mean_dist_yx[ring_rows, ring_cols] * (
                    out_yx[ring_rows, ring_cols] - to_ring
                ) / out_this_new
                new_weighted = j_ring * to_ring / out_this_new
                now_weighted[out_this_new == 0] = 0
                new_weighted[out_this_new == 0] = 0
                mean_dist_yx[ring_rows, ring_cols] = now_weighted + new_weighted
    elif total_unmet_in < 0:
        avail_space = out_yx[ring_rows, ring_cols]
        total_avail = _sum_f(avail_space * (avail_space > 0))
        if total_unmet_in <= -total_avail:
            out_yx[ring_rows, ring_cols] = out_yx[ring_rows, ring_cols] - avail_space * (
                avail_space > 0
            )
            if displaced_yx is not None:
                displaced_yx[ring_rows, ring_cols] = displaced_yx[ring_rows, ring_cols] - avail_space * (
                    avail_space > 0
                )
            total_unmet_out = total_unmet_in + total_avail
        else:
            out_yx[ring_rows, ring_cols] = out_yx[ring_rows, ring_cols] + total_unmet_in * avail_space / total_avail * (
                avail_space > 0
            )
            if displaced_yx is not None:
                displaced_yx[ring_rows, ring_cols] = displaced_yx[ring_rows, ring_cols] + total_unmet_in * avail_space / total_avail * (
                    avail_space > 0
                )
            total_unmet_out = 0

    out_ntrl_yx = in_ntrl_yx - (out_yx - in_yx)
    return out_yx, out_ntrl_yx, total_unmet_out, displaced_yx, mean_dist_yx


def _do_rings_mgmt(
    in_yx: np.ndarray,
    total_unmet_in: float,
    this_area_yx: np.ndarray,
    max_mgmt: float,
    ring_rc: Tuple[np.ndarray, np.ndarray],
    debug: bool = False,
    debug_prefix: str = "",
) -> Tuple[np.ndarray, float]:
    ring_rows, ring_cols = ring_rc
    in_total_ring = in_yx[ring_rows, ring_cols] * this_area_yx[ring_rows, ring_cols]
    out_yx = in_yx.copy()
    total_unmet_out = total_unmet_in

    if total_unmet_in > 0:
        max_tot_ring = max_mgmt * this_area_yx[ring_rows, ring_cols]
        avail = max_tot_ring - in_total_ring
        avail[avail < 0] = 0
        total_avail = _sum_f(avail)
        if debug:
            print(f"{debug_prefix}do_rings: total_unmet_in={total_unmet_in:0.4e}")
            print(f"{debug_prefix}do_rings: total_avail={total_avail:0.4e}")
        if total_avail > 0:
            if total_unmet_in >= total_avail:
                to_ring = avail
                total_unmet_out = total_unmet_in - total_avail
            else:
                to_ring = total_unmet_in * avail / total_avail
                total_unmet_out = 0
            out_total_ring = in_total_ring + to_ring
            if debug:
                print(
                    f"{debug_prefix}do_rings: out_total_ring_sum={_sum_f(out_total_ring):0.4e} "
                    f"unmet_out={total_unmet_out:0.4e}"
                )
            out_ring = out_total_ring / this_area_yx[ring_rows, ring_cols]
            out_ring[this_area_yx[ring_rows, ring_cols] == 0] = 0
            out_yx[ring_rows, ring_cols] = out_ring
    elif total_unmet_in < 0:
        avail = in_total_ring
        total_avail = _sum_f(avail * (avail > 0))
        if debug:
            print(f"{debug_prefix}do_rings: total_unmet_in={total_unmet_in:0.4e}")
            print(f"{debug_prefix}do_rings: total_avail={total_avail:0.4e}")
        if total_unmet_in <= -total_avail:
            out_total_ring = in_total_ring - avail * (avail > 0)
            total_unmet_out = total_unmet_in + total_avail
        else:
            out_total_ring = in_total_ring + total_unmet_in * avail / total_avail * (avail > 0)
            total_unmet_out = 0
        if debug:
            print(
                f"{debug_prefix}do_rings: out_total_ring_sum={_sum_f(out_total_ring):0.4e} "
                f"unmet_out={total_unmet_out:0.4e}"
            )
        out_ring = out_total_ring / this_area_yx[ring_rows, ring_cols]
        out_ring[this_area_yx[ring_rows, ring_cols] == 0] = 0
        out_yx[ring_rows, ring_cols] = out_ring

    return out_yx, total_unmet_out
