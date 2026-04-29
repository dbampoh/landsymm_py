"""Convert harmonized PLUM outputs to LPJ-GUESS-format input files.

Why this module exists
======================
After ``plumharm`` produces year-by-year harmonized scenario outputs
(typically as .mat files containing 2D maps), LPJ-GUESS needs them in
its own text-table format with one row per (lon, lat, year) and one
column per land-cover class / crop / management variable. This module
performs that conversion, producing the four LPJ-GUESS-ready files
referenced by the ``file_lu``, ``file_lucrop``, ``file_Nfert``, and
``file_irrigintens`` parameters in LPJ-GUESS instruction (.ins) files
for each scenario.

The output of this module is what the LandSyMM ↔ LPJ-GUESS coupling
described in Alexander et al. (2018) and Rabin et al. (2020) feeds
into the *final* LPJ-GUESS run that produces ecosystem-service
indicators over the 21st century.

Pseudocode
==========
1) Load PLUMharm + PLUMharm2LPJG options.
2) Resolve harmDirs/forLPJG_dirs.
3) For each scenario:
   - Load yearly .mat outputs.
   - Stack into time arrays.
   - Convert to table format.
   - Save landcover/cropfractions/nfert/irrig.
"""

from __future__ import annotations

import os
from typing import Any, Dict, List, Sequence

import h5py
import numpy as np
from scipy.io import loadmat

from landsymm.common.lpjg_io import (_is_hdf5_mat, _mat_struct_to_dict, maps2table,
                               save_table)

from .plumharm2lpjg_options import PlumHarm2LPJGConfig


class PlumHarm2LPJGPipeline:
    """Class-based PLUMharm2LPJG pipeline (skeleton)."""

    def __init__(self, cfg: PlumHarm2LPJGConfig) -> None:
        self.cfg = cfg

    def run(self) -> None:
        """Run conversion to forLPJG outputs (TODO)."""
        run_plumharm2lpjg(self.cfg)


def run_plumharm2lpjg(cfg: PlumHarm2LPJGConfig) -> None:
    """Run conversion to forLPJG outputs (MATLAB parity)."""
    if cfg.this_dir and not os.path.isdir(cfg.this_dir):
        raise RuntimeError(f"this_dir not found: {cfg.this_dir}")

    plum_dirs = _normalize_dirs(cfg.plum_dirs, cfg.this_dir)
    harm_dirs = cfg.harm_dirs
    if harm_dirs is None:
        harm_dirs = _get_harm_dirs(plum_dirs, cfg.fruitveg_sugar_2oil, cfg.combine_crops)
    harm_dirs = _normalize_dirs(harm_dirs, cfg.this_dir)
    harm_dirs = _check_dirs(harm_dirs, "r")

    if cfg.forLPJG_dirs is None:
        for_lpjg_dirs = []
        for d in harm_dirs:
            for_lpjg_dirs.append(f"{_remove_trailing_slash(d)}.forLPJG")
    else:
        for_lpjg_dirs = list(cfg.forLPJG_dirs)
    for_lpjg_dirs = _check_dirs(for_lpjg_dirs, "rw")

    if cfg.y1_pre is None:
        y1_pre = cfg.year1
    else:
        y1_pre = cfg.y1_pre

    year_list = list(range(cfg.year1, cfg.yearN + 1, cfg.yStep))
    if y1_pre < cfg.year1:
        year_list_xtra = list(range(y1_pre, cfg.year1))
    else:
        year_list_xtra = []

    cf_kgNha_kgNm2 = 1e-4
    mincropfrac = 10 ** (-cfg.out_prec) if cfg.someofall else 0

    for harm_dir, out_dir in zip(harm_dirs, for_lpjg_dirs):
        print(f"Importing harmonization outputs from {harm_dir}...")
        lu_in = {}
        cf_in = {}
        nf_in = {}
        ir_in = {}
        list2map = None

        for yi, this_year in enumerate(year_list):
            print(f"   Reading {this_year}...")
            year_dir = os.path.join(harm_dir, str(this_year))
            s_lu = _load_out_y1(os.path.join(year_dir, f"LandCoverFract.base{cfg.base_year}.mat"))
            s_cf = _load_out_y1(os.path.join(year_dir, f"CropFract.base{cfg.base_year}.mat"))
            s_nf = _load_out_y1(os.path.join(year_dir, f"Fert.base{cfg.base_year}.mat"))
            s_ir = _load_out_y1(os.path.join(year_dir, f"Irrig.base{cfg.base_year}.mat"))

            if yi == 0:
                mask = ~np.isnan(s_lu["maps_YXv"][:, :, 0])
                list2map = np.flatnonzero(mask.ravel(order="F")) + 1
                lu_in = {"varNames": s_lu["varNames"], "yearList": year_list}
                cf_in = {"yearList": year_list}
                nf_in = {"yearList": year_list}
                ir_in = {"yearList": year_list}
                crop_list_rf = s_cf["varNames"]
                crop_list_ir = [f"{c}i" for c in crop_list_rf]
                cf_in["varNames"] = crop_list_ir
                nf_in["varNames"] = crop_list_ir
                ir_in["varNames"] = crop_list_ir
                lu_in["maps_YXvy"] = np.full((*s_lu["maps_YXv"].shape, len(year_list)), np.nan)
                cf_in["maps_YXvy"] = np.full((*s_cf["maps_YXv"].shape, len(year_list)), np.nan)
                nf_in["maps_YXvy"] = np.full((*s_nf["maps_YXv"].shape, len(year_list)), np.nan)
                ir_in["maps_YXvy"] = np.full((*s_ir["maps_YXv"].shape, len(year_list)), np.nan)

            lu_in["maps_YXvy"][:, :, :, yi] = s_lu["maps_YXv"]
            cf_in["maps_YXvy"][:, :, :, yi] = s_cf["maps_YXv"]
            nf_in["maps_YXvy"][:, :, :, yi] = s_nf["maps_YXv"]
            ir_in["maps_YXvy"][:, :, :, yi] = s_ir["maps_YXv"]

        print("Getting arrays...")
        lu_out, lu_header = maps2table(lu_in, list2map)
        cf_out, cf_header = maps2table(cf_in, list2map)
        nf_out, nf_header = maps2table(nf_in, list2map)
        ir_out, ir_header = maps2table(ir_in, list2map)

        print("Checking arrays...")
        for arr, name in [(lu_out, "lu_out"), (cf_out, "cf_out"), (nf_out, "nf_out"), (ir_out, "ir_out")]:
            if np.isnan(arr).any():
                raise RuntimeError(f"Some value of {name} is NaN.")

        ia, _ = _intersect_indices(nf_header, crop_list_ir)
        if np.min(lu_out[:, 3:]) < 0:
            raise RuntimeError("min(lu_out)<0")
        if np.max(lu_out[:, 3:]) >= 1 + 10 ** (-cfg.out_prec):
            raise RuntimeError("max(lu_out)>1")
        if np.max(np.sum(lu_out[:, 3:], axis=1)) >= 1 + 10 ** (-cfg.out_prec):
            raise RuntimeError("max(sum_lu)>1")
        if np.min(cf_out[:, 3:]) < 0:
            raise RuntimeError("min(cf_out)<0")
        if np.max(cf_out[:, 3:]) >= 1 + 10 ** (-cfg.out_prec):
            raise RuntimeError("max(cf_out)>1")
        if np.max(np.sum(cf_out[:, ia], axis=1)) >= 1 + 10 ** (-cfg.out_prec):
            raise RuntimeError("max(sum_cf)>1")
        if np.min(nf_out[:, ia]) < 0:
            raise RuntimeError("min(nf_out)<0")
        if np.min(ir_out[:, ia]) < 0:
            raise RuntimeError("min(ir_out)<0")
        if np.max(ir_out[:, ia]) >= 1 + 10 ** (-cfg.out_prec):
            raise RuntimeError("max(ir_out)>1")

        print("Additional processing...")
        if mincropfrac > 0:
            i_crop = lu_in["varNames"].index("CROPLAND")
            lu_tmp = np.round(lu_out[:, 3:], cfg.out_prec)
            no_cropland = lu_tmp[:, i_crop] < mincropfrac
            i = 0
            while np.any(no_cropland):
                i += 1
                if i > len(cfg.donation_order):
                    raise RuntimeError("GET CROPLAND FROM SOMEWHERE")
                this_donor = cfg.donation_order[i - 1]
                i_this = lu_in["varNames"].index(this_donor)
                involved = (lu_tmp[:, i_this] >= mincropfrac) & no_cropland
                if np.any(involved):
                    transfer_amt = mincropfrac - lu_tmp[involved, i_crop]
                    lu_tmp[involved, i_this] -= transfer_amt
                    lu_tmp[involved, i_crop] += transfer_amt
                    no_cropland = lu_tmp[:, i_crop] == 0
            lu_out[:, 3:] = lu_tmp

            cf_tmp = cf_out[:, ia]
            n_crops = len(ia)
            zero_sum = np.sum(cf_tmp, axis=1) == 0
            cf_tmp[zero_sum, :] = 1 / n_crops
            if np.max(np.sum(cf_tmp, axis=1)) >= 1 + 2 * 10 ** (-cfg.out_prec):
                raise RuntimeError("max(round(sum_cf))>1")
            for c in range(n_crops):
                thiscrop = cf_tmp[:, c]
                is_zero = thiscrop < mincropfrac
                is_max = cf_tmp == np.max(cf_tmp, axis=1)[:, None]
                for i in range(n_crops - 1, 0, -1):
                    tmp = is_max[:, i].copy()
                    sum_left = np.sum(is_max[:, :i], axis=1)
                    tmp[sum_left > 0] = False
                    is_max[:, i] = tmp
                cf_tmp[(is_max & is_zero[:, None])] -= mincropfrac
                cf_tmp[is_zero, c] += mincropfrac
            cf_out[:, ia] = cf_tmp

        if np.min(cf_out[:, 3:]) < mincropfrac:
            raise RuntimeError("min(cf_out)<mincropfrac")
        if np.max(cf_out[:, 3:]) >= 1 + 10 ** (-cfg.out_prec):
            raise RuntimeError("max(cf_out)>1")
        if np.max(np.sum(cf_out[:, 3:], axis=1)) >= 1 + 10 ** (-cfg.out_prec):
            raise RuntimeError("max(sum_cf)>1")
        if np.min(nf_out[:, ia]) < 0:
            raise RuntimeError("min(nf_out)<0")
        if np.min(ir_out[:, ia]) < 0:
            raise RuntimeError("min(ir_out)<0")
        if np.max(ir_out[:, ia]) >= 1 + 10 ** (-cfg.out_prec):
            raise RuntimeError("max(ir_out)>1")

        nf_out[:, ia] = cf_kgNha_kgNm2 * nf_out[:, ia]

        cf_out = np.hstack([cf_out, np.zeros_like(cf_out[:, ia])])
        nf_out = np.hstack([nf_out, np.zeros_like(nf_out[:, ia])])
        ir_out = np.hstack([ir_out, np.zeros_like(ir_out[:, ia])])
        cf_header = cf_header + crop_list_rf
        nf_header = nf_header + crop_list_rf
        ir_header = ir_header + crop_list_rf

        if "ExtraCropi" in cf_header:
            i_extra = cf_header.index("ExtraCrop")
            i_extrai = cf_header.index("ExtraCropi")
            cf_out[:, i_extra] = cf_out[:, i_extrai]
            cf_out = np.delete(cf_out, i_extrai, axis=1)
            nf_out = np.delete(nf_out, i_extrai, axis=1)
            ir_out = np.delete(ir_out, i_extrai, axis=1)
            cf_header.pop(i_extrai)
            nf_header.pop(i_extrai)
            ir_header.pop(i_extrai)
            crop_list_ir = [c for c in crop_list_ir if c != "ExtraCropi"]

        if len(year_list_xtra) > 0:
            lu_out = _prepend_years(lu_out, year_list, year_list_xtra)
            cf_out = _prepend_years(cf_out, year_list, year_list_xtra)
            nf_out = _prepend_years(nf_out, year_list, year_list_xtra)
            ir_out = _prepend_years(ir_out, year_list, year_list_xtra)

        file_out_lu = os.path.join(out_dir, "landcover.txt")
        file_out_cf = os.path.join(out_dir, "cropfractions.txt")
        file_out_nf = os.path.join(out_dir, "nfert.txt")
        file_out_ir = os.path.join(out_dir, "irrig.txt")
        if mincropfrac == 0:
            file_out_lu = file_out_lu.replace(".txt", ".noMinCropFrac.txt")
            file_out_cf = file_out_cf.replace(".txt", ".noMinCropFrac.txt")
            file_out_nf = file_out_nf.replace(".txt", ".noMinCropFrac.txt")
            file_out_ir = file_out_ir.replace(".txt", ".noMinCropFrac.txt")

        print("Saving land cover...")
        save_table(
            lu_header,
            lu_out.astype(np.float32),
            file_out_lu,
            outPrec=cfg.out_prec,
            outWidth=cfg.out_width,
            delimiter=cfg.delimiter,
            overwrite=cfg.overwrite,
            fancy=cfg.fancy,
            save_every_pct=cfg.save_every_pct,
            verbose=cfg.verbose_write,
            gzip_output=cfg.do_gzip,
        )
        print("Saving crop fractions...")
        save_table(
            cf_header,
            cf_out.astype(np.float32),
            file_out_cf,
            outPrec=cfg.out_prec,
            outWidth=cfg.out_width,
            delimiter=cfg.delimiter,
            overwrite=cfg.overwrite,
            fancy=cfg.fancy,
            save_every_pct=cfg.save_every_pct,
            verbose=cfg.verbose_write,
            gzip_output=cfg.do_gzip,
        )

        if mincropfrac == 0:
            print("Not saving fertilization/irrigation because mincropfrac==0.")
        else:
            print("Saving fertilization...")
            save_table(
                nf_header,
                nf_out.astype(np.float32),
                file_out_nf,
                outPrec=cfg.out_prec,
                outWidth=cfg.out_width,
                delimiter=cfg.delimiter,
                overwrite=cfg.overwrite,
                fancy=cfg.fancy,
                save_every_pct=cfg.save_every_pct,
                verbose=cfg.verbose_write,
                gzip_output=cfg.do_gzip,
            )
            print("Saving irrigation...")
            save_table(
                ir_header,
                ir_out.astype(np.float32),
                file_out_ir,
                outPrec=cfg.out_prec,
                outWidth=cfg.out_width,
                delimiter=cfg.delimiter,
                overwrite=cfg.overwrite,
                fancy=cfg.fancy,
                save_every_pct=cfg.save_every_pct,
                verbose=cfg.verbose_write,
                gzip_output=cfg.do_gzip,
            )

    print("All done!")


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


def _intersect_indices(a: List[str], b: List[str]) -> tuple[List[int], List[int]]:
    ia = []
    for idx, val in enumerate(a):
        if val in b:
            ia.append(idx)
    return ia, [b.index(a[i]) for i in ia]


def _prepend_years(arr: np.ndarray, year_list: Sequence[int], year_list_xtra: Sequence[int]) -> np.ndarray:
    n_years = len(year_list)
    n_cells = int(arr.shape[0] / n_years)
    lons = arr[0::n_years, 0]
    lats = arr[0::n_years, 1]
    lons_out = np.repeat(lons, n_years + len(year_list_xtra))
    lats_out = np.repeat(lats, n_years + len(year_list_xtra))
    years_out = np.tile(np.array(list(year_list_xtra) + list(year_list)), n_cells)
    n_vars = arr.shape[1] - 3
    tmp = arr[:, 3:].reshape((n_cells, n_years, n_vars))
    tmp = np.concatenate([np.repeat(tmp[:, :1, :], len(year_list_xtra), axis=1), tmp], axis=1)
    tmp = tmp.reshape((n_cells * (n_years + len(year_list_xtra)), n_vars))
    return np.column_stack([lons_out, lats_out, years_out, tmp])
