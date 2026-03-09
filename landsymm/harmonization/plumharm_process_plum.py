"""Process and reconcile PLUM inputs (MATLAB parity).

Port of PLUMharm_processPLUMin_areaCrops.m.
"""
from __future__ import annotations

import os
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import numpy as np

from landsymm.common.aggregation import aggregate
from landsymm.common.interpolation import inpaint_nans
from landsymm.common.lpjg_io import read_table_then2map


def process_plum_inputs(
    file_in_lcf: str,
    land_area_yx: np.ndarray,
    land_area_2deg_yx: np.ndarray,
    lu_names: Sequence[str],
    bare_frac_y0_yx: np.ndarray | None,
    latest_plum_nfert_2deg_yxv: np.ndarray | None,
    latest_plum_irrig_2deg_yxv: np.ndarray | None,
    plum_to_lpjg: Sequence[str] | None,
    lpjg_crops: Sequence[str],
    norm2extra: float,
    inpaint_method: int | None,
    fruitveg_sugar_2oil: bool,
    allow_unveg: bool,
    out_prec: int = 6,
) -> Tuple[
    Dict[str, Any],
    Dict[str, Any] | None,
    Dict[str, Any] | None,
    Dict[str, Any],
    Dict[str, Any] | None,
    Dict[str, Any] | None,
    np.ndarray | None,
    np.ndarray | None,
    np.ndarray | None,
]:
    """Process PLUM inputs and reconcile crops to LPJ-GUESS."""
    put_urban_here = "BARREN"
    cf_kgNha_kgNm2 = 1e-4
    combine_crops = not plum_to_lpjg
    use_latest_plum_mgmt = latest_plum_nfert_2deg_yxv is not None
    do_interp = inpaint_method is not None

    n_lu = len(lu_names)
    land_area_yxv = np.repeat(land_area_yx[:, :, None], n_lu, axis=2)
    not_bare = np.array([v != "BARREN" for v in lu_names])
    is_crop = np.array([v not in {"NATURAL", "PASTURE", "BARREN"} for v in lu_names])

    # Import LandCoverFract.txt
    s_lcf = read_table_then2map(
        file_in_lcf, verboseIfNoMat=False, force_mat_nosave=True, force_mat_save=False
    )
    s_lcf["varNames"] = _coerce_varnames(s_lcf["varNames"])

    # Move URBAN into BARREN and remove URBAN
    if "URBAN" in s_lcf["varNames"]:
        idx_urban = s_lcf["varNames"].index("URBAN")
        idx_bare = s_lcf["varNames"].index(put_urban_here)
        s_lcf["maps_YXv"][:, :, idx_bare] = np.sum(
            s_lcf["maps_YXv"][:, :, [idx_bare, idx_urban]], axis=2
        )
        s_lcf["maps_YXv"] = np.delete(s_lcf["maps_YXv"], idx_urban, axis=2)
        s_lcf["varNames"].pop(idx_urban)

    # Import detailed LandUse.txt if present
    file_in_dtl = file_in_lcf.replace("LandCoverFract", "LandUse")
    if _file_or_gz_exists(file_in_dtl):
        s_dtl = read_table_then2map(
            file_in_dtl, verboseIfNoMat=False, force_mat_nosave=True, force_mat_save=False
        )
        s_dtl["varNames"] = _coerce_varnames(s_dtl["varNames"])
        is_actual_plum_crops = _contains_all(
            s_dtl["varNames"], "_A", exclude=("ruminants", "monogastrics", "pasture")
        )
        if not combine_crops:
            inds_plum_crops_nfert = _contains_all(
                s_dtl["varNames"], "_FQ", exclude=("ruminants", "monogastrics", "pasture")
            )
            inds_plum_crops_irrig = _contains_all(
                s_dtl["varNames"], "_II", exclude=("ruminants", "monogastrics", "pasture")
            )
        plum_crops = [
            v.replace("_A", "") for v, keep in zip(s_dtl["varNames"], is_actual_plum_crops) if keep
        ]
        s_cropa = {
            "maps_YXv": s_dtl["maps_YXv"][:, :, is_actual_plum_crops]
            * s_lcf["maps_YXv"][:, :, s_lcf["varNames"].index("CROPLAND")][:, :, None]
        }
        if not combine_crops:
            s_nfert = {
                "maps_YXv": _align_mgmt_suffix(
                    s_dtl["maps_YXv"], s_dtl["varNames"], plum_crops, "_FQ"
                )
            }
            s_irrig = {
                "maps_YXv": _align_mgmt_suffix(
                    s_dtl["maps_YXv"], s_dtl["varNames"], plum_crops, "_II"
                )
            }
            if "setaside" not in plum_crops:
                raise RuntimeError(f"setaside not present {file_in_dtl}")
            i_set = plum_crops.index("setaside")
            tmp_n = s_nfert["maps_YXv"][:, :, i_set]
            tmp_i = s_irrig["maps_YXv"][:, :, i_set]
            if np.any(tmp_n > 0) or np.any(tmp_i > 0):
                tmp_n[tmp_n > 0] = 0
                tmp_i[tmp_i > 0] = 0
                s_nfert["maps_YXv"][:, :, i_set] = tmp_n
                s_irrig["maps_YXv"][:, :, i_set] = tmp_i
    else:
        file_in_cropfrac = file_in_lcf.replace("LandCoverFract", "CropFract")
        s_cropf = read_table_then2map(
            file_in_cropfrac, verboseIfNoMat=False, force_mat_nosave=True, force_mat_save=False
        )
        s_cropf["varNames"] = _coerce_varnames(s_cropf["varNames"])
        plum_crops = list(s_cropf["varNames"])
        s_cropa = {
            "maps_YXv": s_cropf["maps_YXv"]
            * s_lcf["maps_YXv"][:, :, s_lcf["varNames"].index("CROPLAND")][:, :, None]
        }
        if not combine_crops:
            file_in_nfert = file_in_lcf.replace("LandCoverFract", "Fert")
            s_nfert = read_table_then2map(
                file_in_nfert, verboseIfNoMat=False, force_mat_nosave=True, force_mat_save=False
            )
            file_in_irrig = file_in_lcf.replace("LandCoverFract", "Irrig")
            s_irrig = read_table_then2map(
                file_in_irrig, verboseIfNoMat=False, force_mat_nosave=True, force_mat_save=False
            )
            s_nfert["varNames"] = _coerce_varnames(s_nfert["varNames"])
            s_irrig["varNames"] = _coerce_varnames(s_irrig["varNames"])
            s_nfert["maps_YXv"] = _align_mgmt_names(
                s_nfert["maps_YXv"], s_nfert["varNames"], plum_crops
            )
            s_irrig["maps_YXv"] = _align_mgmt_names(
                s_irrig["maps_YXv"], s_irrig["varNames"], plum_crops
            )
            if "setaside" not in plum_crops:
                raise RuntimeError(
                    "setaside not present, probably because LandUse.txt wasn't found: "
                    f"{file_in_cropfrac}"
                )

    if combine_crops:
        s_cropa["maps_YXv"] = np.sum(s_cropa["maps_YXv"], axis=2, keepdims=True)
        plum_crops = list(lpjg_crops)
        s_nfert = None
        s_irrig = None
        s_nfert_2deg = None
        s_irrig_2deg = None
        latest_plum_nfert_2deg_yxv = None
        latest_plum_irrig_2deg_yxv = None
        max_orig_nfert = None
    else:
        # Move fruitveg/sugar into Oilcrops if needed
        if fruitveg_sugar_2oil:
            fake_plum = ["fruitveg", "sugar"]
            fake_lpjg = ["Oilcrops", "Oilcrops"]
            for plum_fake, lpjg_real in zip(fake_plum, fake_lpjg):
                if plum_fake not in plum_crops:
                    continue
                if lpjg_real not in lpjg_crops:
                    raise RuntimeError(f"{lpjg_real} not in LPJGcrops")
                plum_real = plum_to_lpjg[lpjg_crops.index(lpjg_real)]
                if plum_real not in plum_crops:
                    raise RuntimeError(f"{plum_real} not in PLUM crops")
                i_real = plum_crops.index(plum_real)
                i_fake = plum_crops.index(plum_fake)
                data_real = s_cropa["maps_YXv"][:, :, i_real]
                data_fake = s_cropa["maps_YXv"][:, :, i_fake]
                wt_real = data_real / (data_real + data_fake)
                wt_fake = data_fake / (data_real + data_fake)
                wt_real[(data_real + data_fake) == 0] = 0
                wt_fake[(data_real + data_fake) == 0] = 0
                s_cropa["maps_YXv"][:, :, i_real] = data_real + data_fake
                s_cropa["maps_YXv"] = np.delete(s_cropa["maps_YXv"], i_fake, axis=2)
                if s_nfert is not None:
                    data_real = s_nfert["maps_YXv"][:, :, i_real]
                    data_fake = s_nfert["maps_YXv"][:, :, i_fake]
                    s_nfert["maps_YXv"][:, :, i_real] = data_real * wt_real + data_fake * wt_fake
                    s_nfert["maps_YXv"] = np.delete(s_nfert["maps_YXv"], i_fake, axis=2)
                if s_irrig is not None:
                    data_real = s_irrig["maps_YXv"][:, :, i_real]
                    data_fake = s_irrig["maps_YXv"][:, :, i_fake]
                    s_irrig["maps_YXv"][:, :, i_real] = data_real * wt_real + data_fake * wt_fake
                    s_irrig["maps_YXv"] = np.delete(s_irrig["maps_YXv"], i_fake, axis=2)
                plum_crops.pop(i_fake)

        # Check crop list compatibility
        plum_to_lpjg_vals = [v for v in plum_to_lpjg if v]
        if plum_crops != [c for c in plum_crops if c in plum_to_lpjg_vals]:
            raise RuntimeError("Mismatch between crops in PLUMcrops and PLUMtoLPJG")

        # Check for fertilization/irrigation on setaside
        if "setaside" in plum_crops:
            i_set = plum_crops.index("setaside")
            tmp_n = s_nfert["maps_YXv"][:, :, i_set]
            tmp_i = s_irrig["maps_YXv"][:, :, i_set]
            if np.any(tmp_n > 0) or np.any(tmp_i > 0):
                raise RuntimeError("Fertilization and/or irrigation on SetAside!")

        # Translate PLUM crop names to LPJ crop names if needed
        if list(lpjg_crops) != list(plum_crops):
            plum_crops = [
                lpjg_crops[plum_to_lpjg.index(crop)] if crop in plum_to_lpjg else crop
                for crop in plum_crops
            ]

        # Move portion of crops into ExtraCrop
        if "ExtraCrop" in plum_crops:
            i_extra = plum_crops.index("ExtraCrop")
            mask_other = np.array([i != i_extra for i in range(len(plum_crops))])
            s_cropa["maps_YXv"][:, :, i_extra] = (
                s_cropa["maps_YXv"][:, :, i_extra]
                + norm2extra * np.sum(s_cropa["maps_YXv"][:, :, mask_other], axis=2)
            )
            s_cropa["maps_YXv"][:, :, mask_other] = s_cropa["maps_YXv"][:, :, mask_other] * (
                1 - norm2extra
            )

    # Concatenate crops and other LUs
    s = {
        "varNames": list(plum_crops)
        + [v for v in s_lcf["varNames"] if v != "CROPLAND"],
        "maps_YXv": np.concatenate(
            [
                s_cropa["maps_YXv"],
                s_lcf["maps_YXv"][:, :, [v != "CROPLAND" for v in s_lcf["varNames"]]],
            ],
            axis=2,
        ),
    }

    # Ensure LU ordering matches baseline
    if list(lu_names) != list(s["varNames"]):
        missing = [v for v in lu_names if v not in s["varNames"]]
        if missing:
            if fruitveg_sugar_2oil and set(missing) == {"FruitAndVeg", "Sugar"}:
                raise RuntimeError(
                    "FruitAndVeg and Sugar are missing from PLUM types; set fruitveg_sugar_2oil to false?"
                )
            raise RuntimeError(f"Types in baseline LU but missing from PLUM: {missing}")
        order = [s["varNames"].index(v) for v in lu_names]
        s["maps_YXv"] = s["maps_YXv"][:, :, order]
        s["varNames"] = list(lu_names)

    # Reorder management to LPJ crop order
    if not combine_crops and list(lpjg_crops) != list(plum_crops):
        order = [plum_crops.index(v) for v in lpjg_crops]
        s_nfert["maps_YXv"] = s_nfert["maps_YXv"][:, :, order]
        s_irrig["maps_YXv"] = s_irrig["maps_YXv"][:, :, order]
    if not combine_crops:
        s_nfert["varNames"] = list(lpjg_crops)
        s_irrig["varNames"] = list(lpjg_crops)

    if not combine_crops:
        s_crop_area = s["maps_YXv"][:, :, is_crop]
        s_nfert["maps_YXv"][s_crop_area == 0] = np.nan
        s_irrig["maps_YXv"][s_crop_area == 0] = np.nan
        s_nfert["maps_YXv"] = s_nfert["maps_YXv"] * cf_kgNha_kgNm2

    # Match mask with overall mask
    s["maps_YXv"][land_area_yxv == 0] = np.nan

    # Harmonize veg/bare fractions
    if bare_frac_y0_yx is not None and np.size(bare_frac_y0_yx) > 0:
        s = _harmonize_vegdbare(s, not_bare, bare_frac_y0_yx, land_area_yx, allow_unveg)

    # Convert to area (m2)
    s["maps_YXv"][np.isnan(s["maps_YXv"])] = 0
    s["maps_YXv"] = s["maps_YXv"] * land_area_yxv

    # Check for bad values
    if not combine_crops:
        tmp_n = s_nfert["maps_YXv"].copy()
        tmp_n[np.isnan(tmp_n)] = 0
        tmp_i = s_irrig["maps_YXv"].copy()
        tmp_i[np.isnan(tmp_i)] = 0
        _check_bad_vals(None, tmp_n, tmp_i, None, lu_names, "import_mgmt", out_prec)

    # Get maps at 2-degree
    if s["maps_YXv"].shape[:2] == (90, 180):
        s_2deg = {"maps_YXv": s["maps_YXv"], "varNames": s["varNames"]}
        s["maps_YXv"] = np.full((360, 720, n_lu), np.nan)
        if not combine_crops:
            s_nfert_2deg = {"maps_YXv": s_nfert["maps_YXv"], "varNames": s_nfert["varNames"]}
            s_irrig_2deg = {"maps_YXv": s_irrig["maps_YXv"], "varNames": s_irrig["varNames"]}
    else:
        s_2deg = {"maps_YXv": aggregate(s["maps_YXv"], 0.5, 2), "varNames": s["varNames"]}
        if not combine_crops:
            s_nfert_2deg = {
                "maps_YXv": _aggregate_mgmt(
                    s_nfert["maps_YXv"], s["maps_YXv"][:, :, is_crop], 0.5, 2
                ),
                "varNames": s_nfert["varNames"],
            }
            s_irrig_2deg = {
                "maps_YXv": _aggregate_mgmt(
                    s_irrig["maps_YXv"], s["maps_YXv"][:, :, is_crop], 0.5, 2
                ),
                "varNames": s_irrig["varNames"],
            }

    if not combine_crops:
        if use_latest_plum_mgmt:
            latest_plum_nfert_2deg_yxv = np.nanmax(
                np.stack([latest_plum_nfert_2deg_yxv, s_nfert_2deg["maps_YXv"]]), axis=0
            )
            latest_plum_irrig_2deg_yxv = np.nanmax(
                np.stack([latest_plum_irrig_2deg_yxv, s_irrig_2deg["maps_YXv"]]), axis=0
            )

        tmp_n = s_nfert_2deg["maps_YXv"].copy()
        tmp_n[np.isnan(tmp_n)] = 0
        tmp_i = s_irrig_2deg["maps_YXv"].copy()
        tmp_i[np.isnan(tmp_i)] = 0
        _check_bad_vals(None, tmp_n, tmp_i, None, lu_names, "import_mgmt", out_prec)

        max_orig_nfert = np.nanmax(s_nfert_2deg["maps_YXv"], axis=(0, 1))
        if do_interp:
            s_nfert_2deg2 = {
                "maps_YXv": _interpolate_mgmt(
                    s_nfert_2deg["maps_YXv"],
                    s_2deg["maps_YXv"][:, :, is_crop],
                    land_area_2deg_yx,
                    lpjg_crops,
                    inpaint_method,
                )
            }
            s_irrig_2deg2 = {
                "maps_YXv": _interpolate_mgmt(
                    s_irrig_2deg["maps_YXv"],
                    s_2deg["maps_YXv"][:, :, is_crop],
                    land_area_2deg_yx,
                    lpjg_crops,
                    inpaint_method,
                )
            }
            s_irrig_2deg2["maps_YXv"][s_irrig_2deg2["maps_YXv"] > 1] = 1
            _check_bad_vals(
                None,
                s_nfert_2deg2["maps_YXv"],
                s_irrig_2deg2["maps_YXv"],
                None,
                lu_names,
                "import_mgmt_2deg2",
                out_prec,
            )
        elif use_latest_plum_mgmt:
            s_nfert_2deg2 = s_nfert_2deg
            s_irrig_2deg2 = s_irrig_2deg

        if use_latest_plum_mgmt:
            repl = (s_2deg["maps_YXv"][:, :, is_crop] == 0) & (
                latest_plum_nfert_2deg_yxv >= 0
            )
            s_nfert_2deg2["maps_YXv"][repl] = latest_plum_nfert_2deg_yxv[repl]
            repl = (s_2deg["maps_YXv"][:, :, is_crop] == 0) & (
                latest_plum_irrig_2deg_yxv >= 0
            )
            s_irrig_2deg2["maps_YXv"][repl] = latest_plum_irrig_2deg_yxv[repl]
            _check_bad_vals(
                None,
                s_nfert_2deg2["maps_YXv"],
                s_irrig_2deg2["maps_YXv"],
                None,
                lu_names,
                "import_mgmt_2deg2.2",
                out_prec,
            )

        if do_interp or use_latest_plum_mgmt:
            s_nfert_2deg["maps_YXv"] = s_nfert_2deg2["maps_YXv"]
            s_irrig_2deg["maps_YXv"] = s_irrig_2deg2["maps_YXv"]
    else:
        s_nfert_2deg = None
        s_irrig_2deg = None
        max_orig_nfert = None

    return (
        s,
        None if combine_crops else s_nfert,
        None if combine_crops else s_irrig,
        s_2deg,
        s_nfert_2deg,
        s_irrig_2deg,
        latest_plum_nfert_2deg_yxv,
        latest_plum_irrig_2deg_yxv,
        max_orig_nfert,
    )


def _file_or_gz_exists(path: str) -> bool:
    return os.path.exists(path) or os.path.exists(f"{path}.gz")


def _coerce_varnames(varnames: Iterable[Any]) -> List[str]:
    if isinstance(varnames, np.ndarray):
        vals = varnames.ravel().tolist()
    else:
        vals = list(varnames)
    out = []
    for v in vals:
        if isinstance(v, bytes):
            out.append(v.decode("utf-8"))
        elif isinstance(v, np.ndarray):
            out.append(str(v.squeeze()))
        else:
            out.append(str(v))
    return out


def _contains_all(
    values: Sequence[str], needle: str, exclude: Sequence[str] = ()
) -> np.ndarray:
    out = []
    for v in values:
        keep = needle in v
        if keep:
            keep = not any(e in v for e in exclude)
        out.append(keep)
    return np.array(out, dtype=bool)


def _align_mgmt_suffix(
    maps_yxv: np.ndarray,
    var_names: Sequence[str],
    plum_crops: Sequence[str],
    suffix: str,
) -> np.ndarray:
    out = np.zeros((maps_yxv.shape[0], maps_yxv.shape[1], len(plum_crops)), dtype=maps_yxv.dtype)
    base_map = {v.replace(suffix, ""): i for i, v in enumerate(var_names) if v.endswith(suffix)}
    for i, crop in enumerate(plum_crops):
        if crop in base_map:
            out[:, :, i] = maps_yxv[:, :, base_map[crop]]
    return out


def _align_mgmt_names(
    maps_yxv: np.ndarray, var_names: Sequence[str], plum_crops: Sequence[str]
) -> np.ndarray:
    out = np.zeros((maps_yxv.shape[0], maps_yxv.shape[1], len(plum_crops)), dtype=maps_yxv.dtype)
    name_map = {v: i for i, v in enumerate(var_names)}
    for i, crop in enumerate(plum_crops):
        if crop in name_map:
            out[:, :, i] = maps_yxv[:, :, name_map[crop]]
    return out


def _harmonize_vegdbare(s: Dict[str, Any], not_bare, bare_frac_y0_yx, land_area_yx, allow_unveg):
    if isinstance(bare_frac_y0_yx, np.ma.MaskedArray):
        bare_frac_y0_yx = np.ma.filled(bare_frac_y0_yx, np.nan)
    bare_frac_y1_yx = np.sum(s["maps_YXv"][:, :, ~not_bare], axis=2)
    if not np.allclose(bare_frac_y0_yx, bare_frac_y1_yx, equal_nan=True, atol=0, rtol=0):
        vegd_frac_y1_yx = np.sum(s["maps_YXv"][:, :, not_bare], axis=2)
        if not allow_unveg and np.any((vegd_frac_y1_yx > 0) & (land_area_yx == 0)):
            raise RuntimeError(
                "PLUM has vegetated fraction where baseline LU has either no land or all unvegetated land"
            )
        vegd_frac_y1_yx_rep = np.repeat(vegd_frac_y1_yx[:, :, None], np.sum(not_bare), axis=2)
        vegd_frac_y1_yxv = s["maps_YXv"][:, :, not_bare] / vegd_frac_y1_yx_rep
        vegd_frac_y1_yxv[vegd_frac_y1_yx_rep == 0] = 0
        bare_frac_y0_yx = np.asarray(bare_frac_y0_yx, dtype=float)
        s["maps_YXv"][:, :, not_bare] = vegd_frac_y1_yxv * (1 - bare_frac_y0_yx[:, :, None])
        s["maps_YXv"][:, :, ~not_bare] = bare_frac_y0_yx[:, :, None]
        maxdiff = np.nanmax(np.abs(np.sum(s["maps_YXv"], axis=2) - 1))
        tol = 1e-12
        if maxdiff > tol:
            raise RuntimeError(
                f"Land use fractions don't sum to 1 within tolerance {tol}; max abs. diff {maxdiff}"
            )
    return s


def _check_bad_vals(in_lu_yxv, in_nfert_yxv, in_irrig_yxv, land_area_yx, lu_names, msg, out_prec):
    is_bare = np.array([v == "BARREN" for v in lu_names])
    if in_lu_yxv is not None and np.any(is_bare):
        vegd_yx = np.sum(in_lu_yxv[:, :, ~is_bare], axis=2)

    if in_lu_yxv is not None and land_area_yx is not None:
        lu_frac = in_lu_yxv / np.repeat(land_area_yx[:, :, None], in_lu_yxv.shape[2], axis=2)
        if np.any(np.isnan(in_lu_yxv)):
            raise RuntimeError(f"NaN(s) in {msg} LU maps!")
        if np.min(lu_frac) < 0:
            raise RuntimeError(f"Negative value(s) in {msg} LU maps! Min {np.min(lu_frac):0.1e}")
        if np.max(lu_frac) > 1 + 10 ** (-out_prec):
            raise RuntimeError(
                f"Value(s) >1 in {msg} LU maps! Max overage {np.max(lu_frac) - 1:0.1e}"
            )
        if np.max(np.sum(lu_frac, axis=2)) > 1 + 10 ** (-out_prec):
            raise RuntimeError(
                f"Sum(s) >1 in {msg} LU maps! Max overage {np.max(np.sum(lu_frac, axis=2)) - 1:0.1e}"
            )
        if np.any(is_bare) and np.any((vegd_yx == 0) & (land_area_yx > 0)):
            raise RuntimeError(
                "Zero vegetation in cell(s) with land! Max land area "
                f"{np.max(land_area_yx[vegd_yx == 0]):0.1e}, total {np.sum(land_area_yx[vegd_yx == 0]):0.1e}."
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


def _aggregate_mgmt(in_array, in_area_array, in_res, out_res):
    res_ratio = int(out_res / in_res)
    in_array = in_array.copy()
    in_array[in_area_array == 0] = np.nan
    in_array[np.isnan(in_array)] = 0
    tmp = np.zeros((in_array.shape[0], in_array.shape[1] // res_ratio) + in_array.shape[2:])
    tmp_area = np.zeros_like(tmp)
    for j in range(res_ratio):
        new_area = in_area_array[:, j::res_ratio, ...]
        tmp_wtd = tmp * tmp_area / (tmp_area + new_area)
        new_wtd = in_array[:, j::res_ratio, ...] * new_area / (tmp_area + new_area)
        mask = (tmp_area + new_area) == 0
        tmp_wtd[mask] = 0
        new_wtd[mask] = 0
        tmp = tmp_wtd + new_wtd
        tmp_area = tmp_area + new_area
    out = np.zeros((in_array.shape[0] // res_ratio, in_array.shape[1] // res_ratio) + in_array.shape[2:])
    out_area = np.zeros_like(out)
    for i in range(res_ratio):
        new_area = tmp_area[i::res_ratio, ...]
        out_wtd = out * out_area / (out_area + new_area)
        new_wtd = tmp[i::res_ratio, ...] * new_area / (out_area + new_area)
        mask = (out_area + new_area) == 0
        out_wtd[mask] = 0
        new_wtd[mask] = 0
        out = out_wtd + new_wtd
        out_area = out_area + new_area
    out[out_area == 0] = np.nan
    return out


def _interpolate_mgmt(in_yxv, area_yxv, land_area_yx, crops, inpaint_method):
    out = np.full_like(in_yxv, np.nan)
    for c, crop in enumerate(crops):
        if crop == "ExtraCrop" or not np.any(in_yxv[:, :, c] > 0):
            out[:, :, c] = np.zeros_like(land_area_yx)
            continue
        out[:, :, c] = inpaint_nans(in_yxv[:, :, c], inpaint_method)
    out[np.repeat(land_area_yx[:, :, None] == 0, len(crops), axis=2)] = 0
    if np.any(out < 0):
        out[out < 0] = 0
    return out
