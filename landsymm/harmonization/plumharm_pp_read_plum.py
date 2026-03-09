"""Read and cache processed PLUM inputs (MATLAB parity)."""
from __future__ import annotations

import fnmatch
import os
from typing import Any, Dict, List, Sequence, Tuple

import h5py
import numpy as np
from scipy.io import loadmat, savemat

from landsymm.common.aggregation import aggregate
from landsymm.common.lpjg_io import (_is_hdf5_mat, _mat_struct_to_dict,
                               read_table_then2map)

from .plumharm_process_plum import process_plum_inputs


def pp_read_plum(
    in_dir: str,
    base_year: int,
    year_list: Sequence[int],
    land_area_yx: np.ndarray,
    lu_names: Sequence[str],
    plum_to_lpjg: Sequence[str] | None,
    lpjg_crops: Sequence[str],
    is_2deg: bool,
    bare_frac_y0_yx: np.ndarray | None,
    norm2extra: float,
    inpaint_method: int | None,
    this_ver: str,
    is_orig: bool,
    fruitveg_sugar_2oil: bool,
    allow_unveg: bool,
) -> Tuple[Dict[str, Any], Dict[str, Any] | None, Dict[str, Any] | None]:
    """Read PLUM inputs, optionally using cached processed MAT file."""
    combine_crops = not plum_to_lpjg
    s_nfert_out = None
    s_irrig_out = None

    matfile = f"{in_dir}.processed.{year_list[0]}-{year_list[-1]}.{this_ver}mat"
    print(matfile)

    if os.path.exists(matfile):
        mat_mtime = os.path.getmtime(matfile)
        txt_mtime = _latest_txt_mtime(in_dir, this_ver)
        if txt_mtime is None or mat_mtime > txt_mtime:
            print("   Loading MAT file...")
            loaded = loadmat(matfile)
            s_out = _mat_to_struct(loaded, "S_out")
            if not combine_crops:
                s_nfert_out = _mat_to_struct(loaded, "S_nfert_out")
                s_irrig_out = _mat_to_struct(loaded, "S_irrig_out")
            return s_out, s_nfert_out, s_irrig_out
        print(".mat file exists but is older than latest .txt file. Regenerating.")

    s_out, s_nfert_out, s_irrig_out = _generate_struct(
        in_dir,
        base_year,
        year_list,
        land_area_yx,
        lu_names,
        plum_to_lpjg,
        lpjg_crops,
        is_2deg,
        bare_frac_y0_yx,
        norm2extra,
        inpaint_method,
        is_orig,
        fruitveg_sugar_2oil,
        allow_unveg,
    )

    print("   Saving MAT file...")
    try:
        if combine_crops:
            savemat(matfile, {"S_out": _prepare_mat_struct(s_out)})
        else:
            savemat(
                matfile,
                {
                    "S_out": _prepare_mat_struct(s_out),
                    "S_nfert_out": _prepare_mat_struct(s_nfert_out),
                    "S_irrig_out": _prepare_mat_struct(s_irrig_out),
                },
            )
    except Exception:
        raise

    print("Done.")
    return s_out, s_nfert_out, s_irrig_out


def _generate_struct(
    in_dir: str,
    base_year: int,
    year_list: Sequence[int],
    land_area_yx: np.ndarray,
    lu_names: Sequence[str],
    plum_to_lpjg: Sequence[str] | None,
    lpjg_crops: Sequence[str],
    is_2deg: bool,
    bare_frac_y0_yx: np.ndarray | None,
    norm2extra: float,
    inpaint_method: int | None,
    is_orig: bool,
    fruitveg_sugar_2oil: bool,
    allow_unveg: bool,
) -> Tuple[Dict[str, Any], Dict[str, Any] | None, Dict[str, Any] | None]:
    combine_crops = not plum_to_lpjg
    s_nfert = None
    s_irrig = None

    s_out: Dict[str, Any] = {"varNames": list(lu_names), "yearList": list(year_list)}
    if not combine_crops:
        s_nfert = {"varNames": list(lpjg_crops), "yearList": list(year_list)}
        s_irrig = {"varNames": list(lpjg_crops), "yearList": list(year_list)}

    n_years = len(year_list)
    land_area_2deg_yx = aggregate(land_area_yx, 0.5, 2) if not is_2deg else land_area_yx

    for y, this_year in enumerate(year_list, start=1):
        file_in_lu = os.path.join(in_dir, str(this_year), "LandCoverFract.mat")
        if not is_orig:
            file_in_lu = file_in_lu.replace(".mat", f".base{base_year}.mat")
        mat_exists = os.path.exists(file_in_lu)
        txt_exists = False
        gz_exists = False
        if not mat_exists:
            file_in_lu = file_in_lu.replace(".mat", ".txt")
            txt_exists = os.path.exists(file_in_lu)
            if not txt_exists:
                gz_exists = os.path.exists(f"{file_in_lu}.gz")

        if not (mat_exists or txt_exists or gz_exists):
            raise RuntimeError(
                f"{this_year} ({file_in_lu.replace('.txt.gz','[.mat,.txt(.gz)]')}) does not exist!"
            )

        if y > 1 and y % 5 == 1:
            print(f"\n   {this_year}... ", end="")
        elif y == 1:
            print(f"   {this_year}... ", end="")
        elif y == n_years:
            print(f"{this_year}...\n", end="")
        else:
            print(f"{this_year}... ", end="")

        if mat_exists:
            lu_tmp_in = _load_out_y1(file_in_lu)
            lu_area_tmp = dict(lu_tmp_in)
            if is_2deg:
                raise RuntimeError("You need to pass in half-degree land area for next step...")
            lu_area_tmp["maps_YXv"] = lu_tmp_in["maps_YXv"] * np.repeat(
                land_area_yx[:, :, None], lu_tmp_in["maps_YXv"].shape[2], axis=2
            )
            if is_2deg:
                lu_area_tmp["maps_YXv"] = aggregate(lu_area_tmp["maps_YXv"], 0.5, 2)

            if combine_crops:
                crop_frac_tmp = {
                    "maps_YXv": np.ones((*lu_area_tmp["maps_YXv"].shape[:2], 1)),
                    "varNames": list(lpjg_crops),
                }
            else:
                file_in_cf = file_in_lu.replace("LandCoverFract", "CropFract")
                crop_frac_tmp = _load_out_y1(file_in_cf)
            crop_area_tmp = dict(crop_frac_tmp)
            n_crops = crop_frac_tmp["maps_YXv"].shape[2]
            i_crop = lu_area_tmp["varNames"].index("CROPLAND")
            crop_area_tmp["maps_YXv"] = crop_frac_tmp["maps_YXv"] * lu_area_tmp["maps_YXv"][
                :, :, i_crop
            ][:, :, None]

            plum_in = {
                "varNames": crop_area_tmp["varNames"]
                + [v for v in lu_area_tmp["varNames"] if v != "CROPLAND"],
                "maps_YXv": np.concatenate(
                    [
                        crop_area_tmp["maps_YXv"],
                        lu_area_tmp["maps_YXv"][:, :, [v != "CROPLAND" for v in lu_area_tmp["varNames"]]],
                    ],
                    axis=2,
                ),
            }

            if not combine_crops:
                file_in_nfert = file_in_lu.replace("LandCoverFract", "Fert")
                nfert_in = _load_out_y1(file_in_nfert)
                file_in_irrig = file_in_lu.replace("LandCoverFract", "Irrig")
                irrig_in = _load_out_y1(file_in_irrig)
                nfert_in["maps_YXv"] = nfert_in["maps_YXv"] * 1e-4
        else:
            plum_in, nfert_in, irrig_in, *_ = process_plum_inputs(
                file_in_lu,
                land_area_yx,
                land_area_2deg_yx,
                lu_names,
                bare_frac_y0_yx,
                None,
                None,
                plum_to_lpjg,
                lpjg_crops,
                norm2extra,
                inpaint_method,
                fruitveg_sugar_2oil,
                allow_unveg,
            )

        if y == 1:
            s_out["maps_YXvy"] = np.full(
                (*plum_in["maps_YXv"].shape[:2], len(lu_names), n_years), np.nan, dtype=float
            )
            if not combine_crops:
                s_nfert["maps_YXvy"] = np.full(
                    (*nfert_in["maps_YXv"].shape[:2], len(lpjg_crops), n_years), np.nan, dtype=float
                )
                s_irrig["maps_YXvy"] = np.full(
                    (*irrig_in["maps_YXv"].shape[:2], len(lpjg_crops), n_years), np.nan, dtype=float
                )

        s_out["maps_YXvy"][:, :, :, y - 1] = plum_in["maps_YXv"]
        if not combine_crops:
            s_nfert["maps_YXvy"][:, :, :, y - 1] = nfert_in["maps_YXv"]
            s_irrig["maps_YXvy"][:, :, :, y - 1] = irrig_in["maps_YXv"]

    if "maps_YXvy" not in s_out:
        raise RuntimeError("No files read!")

    return s_out, s_nfert, s_irrig


def _latest_txt_mtime(in_dir: str, this_ver: str) -> float | None:
    pattern = f"*.{this_ver}txt" if this_ver else "*.txt"
    latest = None
    for root, _dirs, files in os.walk(in_dir):
        for name in files:
            if fnmatch.fnmatch(name, pattern) or fnmatch.fnmatch(name, f"{pattern}.gz"):
                path = os.path.join(root, name)
                mtime = os.path.getmtime(path)
                latest = mtime if latest is None else max(latest, mtime)
    return latest


def _mat_to_struct(loaded: Dict[str, Any], key: str) -> Dict[str, Any]:
    if key not in loaded:
        raise RuntimeError(f"MAT file missing {key}")
    out = _mat_struct_to_dict(loaded[key])
    if "varNames" in out:
        out["varNames"] = _coerce_varnames(out["varNames"])
    if "yearList" in out:
        out["yearList"] = np.array(out["yearList"]).ravel().tolist()
    return out


def _prepare_mat_struct(data: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for k, v in data.items():
        if k == "varNames":
            out[k] = np.array(list(v), dtype=object)
        else:
            out[k] = v
    return out


def _coerce_varnames(varnames: Any) -> List[str]:
    if isinstance(varnames, np.ndarray):
        vals = varnames.ravel().tolist()
    elif isinstance(varnames, (list, tuple)):
        vals = list(varnames)
    else:
        vals = [varnames]
    out = []
    for v in vals:
        if isinstance(v, bytes):
            out.append(v.decode("utf-8"))
        elif isinstance(v, np.ndarray):
            out.append(str(v.squeeze()))
        else:
            out.append(str(v))
    return out


def _load_out_y1(path: str) -> Dict[str, Any]:
    if _is_hdf5_mat(path):
        return _read_out_y1_h5(path)
    loaded = loadmat(path)
    if "out_y1" not in loaded:
        raise RuntimeError("MAT file missing out_y1")
    out = _mat_struct_to_dict(loaded["out_y1"])
    out["varNames"] = _coerce_varnames(out.get("varNames", []))
    return out


def _read_out_y1_h5(path: str) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    with h5py.File(path, "r") as f:
        if "out_y1" not in f:
            raise RuntimeError("HDF5 MAT missing out_y1")
        g = f["out_y1"]
        if "maps_YXv" in g:
            maps = np.array(g["maps_YXv"])
            if "varNames" in g:
                vnames = _read_cellstr_h5(f, g["varNames"])
                if maps.ndim == 3 and maps.shape[0] == len(vnames):
                    maps = np.transpose(maps, (2, 1, 0))
                out["varNames"] = vnames
            out["maps_YXv"] = maps
        if "varNames" in g and "varNames" not in out:
            out["varNames"] = _read_cellstr_h5(f, g["varNames"])
    return out


def _read_cellstr_h5(root: h5py.File, ds: h5py.Dataset) -> List[str]:
    vals = []
    for ref in ds[0]:
        if isinstance(ref, h5py.Reference):
            vals.append(_read_char_h5(root, root[ref]))
        else:
            vals.append(str(ref))
    return vals


def _read_char_h5(root: h5py.File, ds: h5py.Dataset) -> str:
    arr = np.array(ds)
    if arr.dtype == np.uint16:
        return bytes(arr.ravel()).decode("utf-16le").strip("\x00")
    if arr.dtype == np.uint8:
        return bytes(arr.ravel()).decode("utf-8", errors="ignore").strip("\x00")
    return str(arr)
