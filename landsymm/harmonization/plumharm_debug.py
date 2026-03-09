"""Debug output helpers for PLUM harmonization (MATLAB parity)."""
from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np


def _format_sci(value: float, prec: int, zero_as_empty: bool = False) -> str:
    if value == 0 or np.isclose(value, 0):
        return "" if zero_as_empty else "0"
    if np.isnan(value):
        return ""
    return f"{value:.{prec}e}"


def _format_float(value: float, prec: int, zero_as_empty: bool = False) -> str:
    if value == 0 or np.isclose(value, 0):
        return "" if zero_as_empty else "0"
    if np.isnan(value):
        return ""
    return f"{value:.{prec}f}"


def _print_table(headers: Sequence[str], rows: Iterable[Sequence[str]]) -> None:
    rows_list = [list(r) for r in rows]
    widths = [
        max(len(str(h)), max((len(str(r[i])) for r in rows_list), default=0))
        for i, h in enumerate(headers)
    ]
    header_line = "  ".join(str(h).ljust(w) for h, w in zip(headers, widths))
    print(header_line)
    for row in rows_list:
        print("  ".join(str(c).ljust(w) for c, w in zip(row, widths)))


def debug_out(
    debug_title: str,
    debug_kind: str,
    agri_yxv: np.ndarray,
    land_area_yx: np.ndarray,
    debug_ij: Sequence[int],
    lu_names: Sequence[str],
    out_prec: int = 2,
) -> None:
    """Mirror PLUMharm_debugOut.m."""
    if debug_ij is None:
        return
    i = debug_ij[0] - 1
    j = debug_ij[1] - 1
    nlu = len(lu_names)

    rows = []
    for v in range(nlu):
        area = float(agri_yxv[i, j, v])
        area_str = _format_sci(area, out_prec)
        if land_area_yx.ndim == 3:
            this_area = float(land_area_yx[i, j, min(v, land_area_yx.shape[2] - 1)])
        else:
            this_area = float(land_area_yx[i, j])
        frac = area / this_area if this_area != 0 else 0.0
        frac_str = _format_sci(frac, out_prec)
        rows.append([lu_names[v], area_str, frac_str])

    if debug_kind == "areas":
        is_cropland = np.array(
            [n not in ("PASTURE", "NATURAL", "BARREN") for n in lu_names], dtype=bool
        )
        is_vegd = np.array([n != "BARREN" for n in lu_names], dtype=bool)
        cropland_area = float(np.sum(agri_yxv[i, j, is_cropland]))
        vegd_area = float(np.sum(agri_yxv[i, j, is_vegd]))
        if land_area_yx.ndim == 3:
            denom = float(np.squeeze(land_area_yx[i, j, 0]))
        else:
            denom = float(land_area_yx[i, j])
        cropland_frac = cropland_area / denom if denom != 0 else 0.0
        vegd_frac = vegd_area / denom if denom != 0 else 0.0
        rows.append(["CROPLAND", _format_sci(cropland_area, out_prec), _format_sci(cropland_frac, out_prec)])
        rows.append(["vegd", _format_sci(vegd_area, out_prec), _format_float(vegd_frac, out_prec)])

    _print_table([debug_title, debug_kind, "Frac"], rows)


def debug_out_deltas(
    debug_title: str,
    debug_kind: str,
    y0_yxv: np.ndarray,
    y1_yxv: np.ndarray,
    debug_ij: Sequence[int],
    lu_names: Sequence[str],
    out_prec: int = 2,
) -> None:
    """Mirror PLUMharm_debugOut_deltas.m."""
    if debug_ij is None:
        return
    i = debug_ij[0] - 1
    j = debug_ij[1] - 1
    nlu = len(lu_names)

    areas_y0 = np.zeros(nlu + 2)
    areas_y1 = np.zeros(nlu + 2)
    areas_y0[:nlu] = y0_yxv[i, j, :]
    areas_y1[:nlu] = y1_yxv[i, j, :]
    is_cropland = np.array([n not in ("PASTURE", "NATURAL", "BARREN") for n in lu_names], dtype=bool)
    is_vegd = np.array([n != "BARREN" for n in lu_names], dtype=bool)
    areas_y0[-2] = float(np.sum(y0_yxv[i, j, is_cropland]))
    areas_y1[-2] = float(np.sum(y1_yxv[i, j, is_cropland]))
    areas_y0[-1] = float(np.sum(y0_yxv[i, j, is_vegd]))
    areas_y1[-1] = float(np.sum(y1_yxv[i, j, is_vegd]))

    diffs = areas_y1 - areas_y0
    with np.errstate(divide="ignore", invalid="ignore"):
        diffs_pct = diffs / areas_y0 * 100

    rows = []
    names = list(lu_names) + ["CROPLAND", "vegd"]
    for name, a0, a1, d, dp in zip(names, areas_y0, areas_y1, diffs, diffs_pct):
        a0_str = _format_sci(float(a0), out_prec)
        a1_str = _format_sci(float(a1), out_prec)
        d_str = _format_sci(float(d), out_prec, zero_as_empty=True)
        if d > 0 and d_str:
            d_str = f"+{d_str}"
        dp_str = _format_float(float(dp), out_prec, zero_as_empty=True)
        if d > 0 and dp_str:
            dp_str = f"+{dp_str}"
        if dp_str:
            dp_str = f"{dp_str} %"
        rows.append([name, a0_str, a1_str, d_str, dp_str])

    _print_table(
        [debug_title, f"{debug_kind}_y0", f"{debug_kind}_y1", "diffs", "diffs_pct"],
        rows,
    )
