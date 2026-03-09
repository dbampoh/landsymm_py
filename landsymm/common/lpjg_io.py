"""LPJ-GUESS I/O helpers (MATLAB-parity focused)."""
# TODO:
# - implement HDF5 (.mat v7.3) cache read/write if needed
from __future__ import annotations

import gzip
import os
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import h5py
import numpy as np
import pandas as pd
from scipy.io import loadmat, savemat


def read_table(
    path: str,
    *,
    dont_save_MAT: bool = False,
    do_save_MAT: bool = True,
    verbose: bool = False,
    verboseIfNoMat: bool = True,
    dispPrefix: str = "",
    force_as_gridlist: bool = False,
) -> pd.DataFrame:
    """Read LPJ-GUESS table with MATLAB-like semantics."""
    if dont_save_MAT and do_save_MAT:
        warnings.warn("dont_save_MAT and do_save_MAT can't both be true. Assuming no save.")
        do_save_MAT = False
    if verboseIfNoMat and verbose:
        verbose = False

    in_file = _resolve_symlink(path)

    # If file doesn't exist, check for gz (including symlink target).
    in_file_gz = f"{in_file}.gz"
    in_file_gz_target = _resolve_symlink(in_file_gz)
    if in_file_gz != in_file_gz_target:
        in_file = in_file_gz_target

    # Expand wildcard
    if "*" in in_file:
        matches = sorted(Path().glob(in_file))
        if len(matches) == 0:
            raise FileNotFoundError(f"No match found for {in_file}")
        if len(matches) > 1:
            raise RuntimeError(f"More than one match found for {in_file}")
        tmp = str(matches[0].resolve())
        if verbose or verboseIfNoMat:
            print(f"Resolving\n{in_file}\ninto\n{tmp}")
        in_file = tmp

    # Strip .gz if provided
    if in_file.endswith(".gz"):
        in_file = in_file[:-3]

    in_matfile = f"{in_file}.mat"
    if os.path.exists(in_matfile):
        if verbose:
            print(f"{dispPrefix}   Loading table from MAT-file...")
        try:
            mat_table = _load_table_from_mat(in_matfile)
            cols = set(str(c) for c in mat_table.columns)
            if {"Lon", "Lat"} & cols:
                return mat_table
            # If MAT table lacks headers, fall back to text read.
            warnings.warn("MAT table missing Lon/Lat headers; reading text instead.")
        except Exception as exc:
            warnings.warn(
                "Problem loading MAT file. Will read text instead.\n"
                f"Intercepted error message follows:\n{exc}"
            )

    did_unzip = _gunzip_if_needed(in_file, verbose or verboseIfNoMat, dispPrefix)
    if verbose or verboseIfNoMat:
        print(f"{dispPrefix}   Making table...")

    # Error if empty
    if os.path.exists(in_file) and os.path.getsize(in_file) == 0:
        raise RuntimeError(f"{in_file}: file is empty!")

    out_table = _import_to_table(in_file, verbose or verboseIfNoMat, dispPrefix, force_as_gridlist)
    if not dont_save_MAT:
        if do_save_MAT:
            try:
                savemat(in_matfile, {"out_table": out_table})
            except Exception:
                warnings.warn(f"Unable to save {in_matfile}")
        else:
            # No interactive prompt; skip save if not forced.
            pass

    if did_unzip:
        if verbose or verboseIfNoMat:
            print(f"{dispPrefix}      Deleting unzipped in_file...")
        os.remove(in_file)

    return _standardize_colnames(out_table)


def read_table_then2map(
    path: str,
    *,
    xres: float = np.nan,
    yres: float = np.nan,
    lat_orient: str = "",
    lon_orient: str = "",
    verbose: bool = False,
    quiet: bool = False,
    verboseIfNoMat: bool = False,
    force_mat_save: bool = True,
    force_mat_nosave: bool = False,
    list_to_map_in: Sequence[int] | None = None,
    dataType: str = "double",
    in_prec: int = 2,
    force_as_gridlist: bool = False,
    drop_northpole: bool = False,
    drop_southpole: bool = False,
    lons_centered_on_180: bool = False,
) -> dict:
    """Read table and convert to map (MATLAB-parity)."""
    if force_mat_nosave and force_mat_save:
        warnings.warn("force_mat_save and force_mat_nosave can't both be true. MATfile will not be saved.")
        force_mat_save = False
    if verboseIfNoMat and verbose:
        verbose = False

    in_file = _strip_known_suffixes(path)
    in_file = _resolve_symlink(in_file)
    in_file = _expand_wildcard(in_file)

    # Optional cache (best-effort)
    in_matfile_maps = f"{in_file}.maps.mat"
    if not force_mat_nosave and os.path.exists(in_matfile_maps):
        try:
            if _is_hdf5_mat(in_matfile_maps):
                cached = _read_maps_mat_h5(in_matfile_maps)
                if cached:
                    return cached
            else:
                loaded = loadmat(in_matfile_maps)
                if "out_struct" in loaded:
                    cached = _mat_struct_to_dict(loaded["out_struct"])
                    list_to_map_cached = np.asarray(cached.get("list_to_map", [])).ravel()
                    if list_to_map_cached.size > 1 or "maps_YXv" in cached or "maps_YXvy" in cached:
                        return cached
        except Exception:
            pass

    if verboseIfNoMat or verbose:
        name = Path(in_file).name
        print(f"{name}:")
        print("   Making maps...")

    in_table = read_table(
        in_file,
        verbose=verbose,
        verboseIfNoMat=verboseIfNoMat,
        dont_save_MAT=force_mat_nosave,
        do_save_MAT=force_mat_save,
        force_as_gridlist=force_as_gridlist,
    )
    is_gridlist = in_table.shape[1] == 2

    yearList, multi_yrs, varNames, list_to_map, xres, yres, found, lat_extent = _get_indices(
        in_table,
        xres,
        yres,
        list_to_map_in,
        lat_orient,
        lon_orient,
        drop_northpole,
        drop_southpole,
        lons_centered_on_180,
        verboseIfNoMat,
        verbose,
        in_prec,
    )
    list_to_map = np.asarray(list_to_map).ravel()

    if is_gridlist:
        mask_YX = _make_maps(
            xres,
            yres,
            multi_yrs,
            yearList,
            varNames,
            in_table,
            list_to_map,
            is_gridlist,
            found,
            lat_extent,
            dataType,
            verboseIfNoMat,
            verbose,
            quiet,
        )
    else:
        if multi_yrs:
            maps_YXvy = _make_maps(
                xres,
                yres,
                multi_yrs,
                yearList,
                varNames,
                in_table,
                list_to_map,
                is_gridlist,
                found,
                lat_extent,
                dataType,
                verboseIfNoMat,
                verbose,
                quiet,
            )
        else:
            maps_YXv = _make_maps(
                xres,
                yres,
                multi_yrs,
                yearList,
                varNames,
                in_table,
                list_to_map,
                is_gridlist,
                found,
                lat_extent,
                dataType,
                verboseIfNoMat,
                verbose,
                quiet,
            )

    # lonlats
    nlatdeg = lat_extent[1] - lat_extent[0]
    nlon = _get_nlonlat_in_tolerance(360, xres, quiet)
    nlat = _get_nlonlat_in_tolerance(nlatdeg, yres, quiet)
    lons_map_yx = np.tile(np.arange(-180 + xres / 2, 180 + 1e-9, xres), (nlat, 1))
    lats_map_yx = np.tile(
        np.arange(lat_extent[0] + yres / 2, lat_extent[1] + 1e-9, yres).reshape(-1, 1),
        (1, nlon),
    )
    lons_out = lons_map_yx.ravel(order="F")[list_to_map - 1]
    lats_out = lats_map_yx.ravel(order="F")[list_to_map - 1]
    lonlats = np.column_stack([lons_out, lats_out])

    out_struct: dict = {
        "list_to_map": list_to_map,
        "lonlats": lonlats,
        "lat_extent": lat_extent,
        "lat_orient": lat_orient,
    }
    if is_gridlist:
        out_struct["mask_YX"] = mask_YX
    else:
        out_struct["varNames"] = list(varNames)
        if multi_yrs:
            out_struct["maps_YXvy"] = maps_YXvy
            out_struct["yearList"] = np.array(yearList)
        else:
            out_struct["maps_YXv"] = maps_YXv

    if not force_mat_nosave:
        try:
            savemat(in_matfile_maps, {"out_struct": out_struct})
        except Exception:
            pass
    if verboseIfNoMat or verbose:
        print("   Done.")
    return out_struct


def read2geoarray(
    path: str,
    *,
    xres: float = np.nan,
    yres: float = np.nan,
    lat_orient: str = "",
    lon_orient: str = "",
    verbose: bool = False,
    verboseIfNoMat: bool = True,
    force_mat_save: bool = True,
    force_mat_nosave: bool = False,
    dataType: str = "double",
    in_prec: int = 2,
    target=None,
    target_lat_orient: str = "",
    target_lon_orient: str = "",
    trimfirstyear_ifneeded: bool = False,
    drop_northpole: bool = False,
    drop_southpole: bool = False,
    lons_centered_on_180: bool = False,
    fill_missing_cells_with_nan: bool = False,
    gridlist_warn: bool = True,
) -> dict:
    """Read to GeoArray-like structure (MATLAB-parity)."""
    if force_mat_nosave and force_mat_save:
        warnings.warn("force_mat_save and force_mat_nosave can't both be true. MATfile will not be saved.")
        force_mat_save = False
    if verboseIfNoMat and verbose:
        verbose = False

    in_file = _strip_known_suffixes(path)
    in_file = _resolve_symlink(in_file)
    in_file = _expand_wildcard(in_file)

    # Try .garr.mat (best-effort)
    in_matfile_garr = f"{in_file}.garr.mat"
    if os.path.exists(in_matfile_garr) and target is None:
        try:
            loaded = loadmat(in_matfile_garr)
            if "out_struct" in loaded:
                out_struct = _mat_struct_to_dict(loaded["out_struct"])
                return _normalize_garr_struct(out_struct)
        except Exception:
            pass

    # MATLAB parity path: when target is provided and table has years, build garr
    if target is not None:
        try:
            table_in = read_table(
                in_file,
                verboseIfNoMat=verboseIfNoMat,
                do_save_MAT=False,
                force_as_gridlist=False,
            )
            if "Year" in table_in.columns and len(np.unique(table_in["Year"].to_numpy())) > 1:
                return _table_to_garr_with_target(table_in, target)
        except Exception:
            pass

    out_struct = read_table_then2map(
        in_file,
        xres=xres,
        yres=yres,
        lat_orient=lat_orient,
        lon_orient=lon_orient,
        verbose=verbose,
        verboseIfNoMat=verboseIfNoMat,
        force_mat_save=force_mat_save,
        force_mat_nosave=force_mat_nosave,
        dataType=dataType,
        in_prec=in_prec,
        drop_northpole=drop_northpole,
        drop_southpole=drop_southpole,
        lons_centered_on_180=lons_centered_on_180,
    )
    # Gridlist case: return mask/list2map/lonlats
    if "mask_YX" in out_struct:
        out_struct["list2map"] = np.asarray(out_struct.pop("list_to_map")).ravel()
        lat_extent = out_struct.get("lat_extent", [-90.0, 90.0])
        lat_extent = np.ravel(lat_extent).tolist()
        if len(lat_extent) != 2:
            raise RuntimeError(f"lat_extent malformed: {lat_extent}")
        nlat = out_struct["mask_YX"].shape[0]
        nlon = out_struct["mask_YX"].shape[1]
        xres2 = 360 / nlon
        yres2 = (lat_extent[1] - lat_extent[0]) / nlat
        lons_map_yx = np.tile(np.arange(-180 + xres2 / 2, 180 + 1e-9, xres2), (nlat, 1))
        lats_map_yx = np.tile(
            np.arange(lat_extent[0] + yres2 / 2, lat_extent[1] + 1e-9, yres2).reshape(-1, 1),
            (1, nlon),
        )
        lons_out = lons_map_yx.ravel(order="F")[out_struct["list2map"] - 1]
        lats_out = lats_map_yx.ravel(order="F")[out_struct["list2map"] - 1]
        out_struct["lonlats"] = np.column_stack([lons_out, lats_out])
        return out_struct

    # If target provided, convert maps to garray using target list2map directly
    if target is not None and ("maps_YXv" in out_struct or "maps_YXvy" in out_struct):
        if isinstance(target, dict):
            lonlats_target = target.get("lonlats")
            list2map_target = target.get("list2map")
        elif isinstance(target, (list, tuple)) and len(target) == 2:
            lonlats_target, list2map_target = target
        else:
            raise RuntimeError("Input target is of invalid class")
        if lonlats_target is None or list2map_target is None:
            raise RuntimeError("target must provide lonlats and list2map")

        lonlats_target = np.asarray(lonlats_target)
        list2map_target = np.asarray(list2map_target, dtype=int).ravel()

        if "maps_YXvy" in out_struct:
            maps = out_struct["maps_YXvy"]
            list2map_all = _get_list2map_all(list2map_target, maps.shape[:2], maps.shape[2:])
            flat = maps.ravel(order="F")[list2map_all - 1]
            garr_xvy = flat.reshape((len(list2map_target), maps.shape[2], maps.shape[3]), order="F")
            out_struct = {
                "list2map": list2map_target,
                "lonlats": lonlats_target,
                "varNames": out_struct.get("varNames", []),
                "yearList": out_struct.get("yearList", []),
                "map_size": maps.shape[:2],
                "garr_xvy": garr_xvy,
            }
        else:
            maps = out_struct["maps_YXv"]
            list2map_all = _get_list2map_all(list2map_target, maps.shape[:2], maps.shape[2:])
            flat = maps.ravel(order="F")[list2map_all - 1]
            garr_xv = flat.reshape((len(list2map_target), maps.shape[2]), order="F")
            out_struct = {
                "list2map": list2map_target,
                "lonlats": lonlats_target,
                "varNames": out_struct.get("varNames", []),
                "map_size": maps.shape[:2],
                "garr_xv": garr_xv,
            }
    else:
        # Convert maps to garray-like representation
        out_struct = _maps_to_garr(out_struct)
    # Ensure lonlats for garr outputs
    if "lonlats" not in out_struct:
        lat_extent = out_struct.get("lat_extent", [-90.0, 90.0])
        map_size = out_struct.get("map_size")
        if map_size is not None:
            nlat, nlon = map_size
            xres2 = 360 / nlon
            yres2 = (lat_extent[1] - lat_extent[0]) / nlat
            lons_map_yx = np.tile(np.arange(-180 + xres2 / 2, 180 + 1e-9, xres2), (nlat, 1))
            lats_map_yx = np.tile(
                np.arange(lat_extent[0] + yres2 / 2, lat_extent[1] + 1e-9, yres2).reshape(-1, 1),
                (1, nlon),
            )
            lons_out = lons_map_yx.ravel(order="F")[out_struct["list2map"] - 1]
            lats_out = lats_map_yx.ravel(order="F")[out_struct["list2map"] - 1]
            out_struct["lonlats"] = np.column_stack([lons_out, lats_out])

    # Align to target gridlist if provided (for non-map garr inputs)
    if target is not None and ("garr_xv" in out_struct or "garr_xvy" in out_struct):
        if isinstance(target, dict):
            lonlats_target = target.get("lonlats")
            list2map_target = target.get("list2map")
        elif isinstance(target, (list, tuple)) and len(target) == 2:
            lonlats_target, list2map_target = target
        else:
            raise RuntimeError("Input target is of invalid class")

        if lonlats_target is None or list2map_target is None:
            raise RuntimeError("target must provide lonlats and list2map")

        lonlats_target = np.asarray(lonlats_target)
        list2map_target = np.asarray(list2map_target, dtype=int)

        if "garr_xvy" in out_struct:
            new_garr = np.full((len(list2map_target), out_struct["garr_xvy"].shape[1], out_struct["garr_xvy"].shape[2]), np.nan)
        elif "garr_xv" in out_struct:
            new_garr = np.full((len(list2map_target), out_struct["garr_xv"].shape[1]), np.nan)
        else:
            raise RuntimeError("in_struct does not appear to include either garr_xvy or garr_xv")

        in_lonlats = out_struct["lonlats"]
        in_index = {tuple(row): idx for idx, row in enumerate(in_lonlats)}
        ia: list[int] = []
        ib: list[int] = []
        for t_idx, row in enumerate(lonlats_target):
            key = tuple(row)
            if key in in_index:
                ia.append(t_idx)
                ib.append(in_index[key])
        ia_idx = np.array(ia, dtype=int)
        ib_idx = np.array(ib, dtype=int)
        if "garr_xvy" in out_struct:
            if ia_idx.size:
                new_garr[ia_idx, :, :] = out_struct["garr_xvy"][ib_idx, :, :]
            out_struct["garr_xvy"] = new_garr
        else:
            if ia_idx.size:
                new_garr[ia_idx, :] = out_struct["garr_xv"][ib_idx, :]
            out_struct["garr_xv"] = new_garr

        out_struct["lonlats"] = lonlats_target
        out_struct["list2map"] = list2map_target

    return out_struct


def maps2table(in_struct: dict, list2map: Sequence[int]) -> tuple[np.ndarray, list[str]]:
    """Convert maps to table (MATLAB-parity)."""
    if "maps_YXvy" in in_struct and "yearList" in in_struct:
        with_years = True
        in_maps = in_struct["maps_YXvy"].astype(float)
        year_list = np.array(in_struct["yearList"])
    elif "maps_YXv" in in_struct and "yearList" not in in_struct:
        with_years = False
        in_maps = in_struct["maps_YXv"].astype(float)
        year_list = None
    else:
        raise RuntimeError("Invalid in_struct: missing maps_YXv(y)")

    nlat, nlon = in_maps.shape[:2]
    xres = 360 / nlon
    yres = 180 / nlat
    lons = np.arange(-180 + xres / 2, 180 - xres / 2 + 1e-9, xres)
    lats = np.arange(-90 + yres / 2, 90 - yres / 2 + 1e-9, yres)
    lons_map = np.tile(lons, (len(lats), 1))
    lats_map = np.tile(lats.reshape(-1, 1), (1, len(lons)))

    list2map = np.array(list2map, dtype=int)
    ncells_out = len(list2map)
    nvars = in_maps.shape[2]

    if with_years:
        nyears = len(year_list)
        out_array = np.full((ncells_out * nyears, 3 + nvars), np.nan, dtype=float)
    else:
        out_array = np.full((ncells_out, 2 + nvars), np.nan, dtype=float)

    for c, idx in enumerate(list2map):
        i, j = np.unravel_index(idx - 1, (nlat, nlon), order="F")
        this_lon = lons_map[i, j]
        this_lat = lats_map[i, j]
        if with_years:
            this_cell_yv = np.transpose(in_maps[i, j, :, :], (1, 0))
            lonlatyr = np.column_stack(
                [np.full(nyears, this_lon), np.full(nyears, this_lat), year_list]
            )
            these_rows = np.column_stack([lonlatyr, this_cell_yv])
            c1 = c * nyears
            cN = c1 + nyears
            out_array[c1:cN, :] = these_rows
        else:
            out_array[c, :] = np.concatenate([[this_lon, this_lat], in_maps[i, j, :]])

    if np.isnan(out_array).any():
        raise RuntimeError("At least one member of out_array is NaN!")

    if with_years:
        out_header_cell = ["Lon", "Lat", "Year"] + list(in_struct["varNames"])
    else:
        out_header_cell = ["Lon", "Lat"] + list(in_struct["varNames"])

    return out_array, out_header_cell


def save_table(
    header,
    out_array,
    out_file: str,
    *,
    outPrec: int = 3,
    outPrec_lonlat: int = 2,
    outWidth: int = 5,
    varNameAlign: str = "L",
    dataAlign: str = "L",
    fancy: bool = False,
    progress_step_pct: int = 5,
    save_every_n: int = 1000,
    save_every_pct: int = 5,
    delimiter: str = " ",
    overwrite: bool = False,
    verbose: bool = True,
    justZeroCols: Sequence[int] | None = None,
    gzip_output: bool = False,
) -> None:
    """Save table with MATLAB-like formatting."""
    if not overwrite and os.path.exists(out_file):
        warnings.warn(f"{out_file} already exists. Specify overwrite=true to overwrite. Skipping.")
        return

    if justZeroCols is None:
        justZeroCols = []

    in_header_cell, in_header_str = _get_header_cells(header, out_array)

    if fancy:
        # Fancy writer is rarely used; fall back to fast writer for parity.
        fancy = False

    _write_fast(
        in_header_str,
        in_header_cell,
        out_file,
        out_array,
        outPrec,
        outPrec_lonlat,
        outWidth,
        delimiter,
        verbose,
        justZeroCols,
        save_every_pct,
    )

    if gzip_output:
        if verbose:
            print("Zipping...")
        with open(out_file, "rb") as f_in, gzip.open(f"{out_file}.gz", "wb") as f_out:
            f_out.writelines(f_in)
        os.remove(out_file)
    if verbose:
        print("Done.")


def _write_fast(
    in_header_str: str | None,
    in_header_cell: list[str],
    out_file: str,
    out_array,
    outPrec: int,
    outPrec_lonlat: int,
    outWidth: int,
    delimiter: str,
    verbose: bool,
    justZeroCols: Sequence[int],
    save_every_pct: int,
) -> None:
    out_header_str = in_header_str or _format_header(in_header_cell, outWidth, delimiter)

    i_lat, i_lon, i_year = _get_lat_lon_year_cols(in_header_cell)

    fmt_cols: list[str] = []
    int_cols: list[bool] = []
    for ii in range(len(in_header_cell)):
        if ii > 0:
            fmt_cols.append(delimiter)
        if ii == i_lat or ii == i_lon:
            fmt_cols.append(f"%-{outWidth}.{outPrec_lonlat}f")
            int_cols.append(False)
        elif ii == i_year or (ii + 1) in justZeroCols:
            fmt_cols.append(f"%-{outWidth}u")
            int_cols.append(True)
        else:
            fmt_cols.append(f"%-{outWidth}.{outPrec}f")
            int_cols.append(False)
    row_fmt = "".join(fmt_cols) + " \n"

    with open(out_file, "w") as f:
        f.write(f"{out_header_str} \n")

    if isinstance(out_array, dict):
        n_cells = len(out_array["list2map"])
        chunk_size = max(int(np.floor(save_every_pct / 100 * n_cells)), 1)
        n_chunks = int(np.ceil(n_cells / chunk_size))

        for ii in range(n_chunks):
            with open(out_file, "a") as f:
                c1 = ii * chunk_size
                cN = min((ii + 1) * chunk_size, n_cells)
                this_chunk = cN - c1

                if "garr_xvy" in out_array:
                    nyears = len(out_array["yearList"])
                    col_lons = np.repeat(out_array["lonlats"][c1:cN, 0], nyears)
                    col_lats = np.repeat(out_array["lonlats"][c1:cN, 1], nyears)
                    col_lonlats = np.column_stack([col_lons, col_lats])
                    col_years = np.tile(np.array(out_array["yearList"]), this_chunk)

                    tmp = out_array["garr_xvy"][c1:cN, :, :]
                    tmp = np.transpose(tmp, (0, 2, 1))
                    tmp = tmp.reshape((this_chunk * nyears, tmp.shape[2]), order="C")
                    out_array_tmp = np.column_stack([col_lonlats, col_years, tmp])
                elif "garr_xv" in out_array:
                    out_array_tmp = np.column_stack([out_array["lonlats"][c1:cN, :], out_array["garr_xv"][c1:cN, :]])
                else:
                    raise RuntimeError("in_data has neither garr_xvy nor garr_xv")

                out_array_tmp[out_array_tmp == 0] = 0
                for row in out_array_tmp:
                    row_vals = [int(v) if int_cols[i] else float(v) for i, v in enumerate(row)]
                    f.write(row_fmt % tuple(row_vals))

            if verbose:
                this_pct = round((ii + 1) / n_chunks * 100)
                print(f"{this_pct}%...")
    else:
        if isinstance(out_array, pd.DataFrame):
            out_array = out_array.to_numpy()
        out_array = np.asarray(out_array).T
        nrows = out_array.shape[1]
        n_chunks = int(np.ceil(100 / save_every_pct))
        chunk_size = int(np.ceil(nrows * 100 / save_every_pct))
        for ii in range(n_chunks):
            with open(out_file, "a") as f:
                c1 = ii * chunk_size
                cN = min((ii + 1) * chunk_size, nrows)
                out_array_tmp = out_array[:, c1:cN]
                out_array_tmp[out_array_tmp == 0] = 0
                for row in out_array_tmp.T:
                    row_vals = [int(v) if int_cols[i] else float(v) for i, v in enumerate(row)]
                    f.write(row_fmt % tuple(row_vals))
            if verbose and ii != n_chunks - 1:
                this_pct = round((ii + 1) / n_chunks * 100)
                print(f"{this_pct}%...")


def _format_header(header_cell: Sequence[str], outWidth: int, delimiter: str) -> str:
    var_name_max = max(len(x) for x in header_cell)
    tmp = []
    for name in header_cell:
        if name in {"Lon", "Lat", "Year"}:
            width = 5
        else:
            width = min(outWidth, var_name_max)
        tmp.append(name.ljust(width))
    return delimiter.join(tmp)


def _get_header_cells(header, out_array) -> tuple[list[str], str | None]:
    if isinstance(out_array, dict):
        if "garr_xvy" in out_array:
            in_header_cell = ["Lon", "Lat", "Year"] + list(out_array["varNames"])
        elif "garr_xv" in out_array:
            in_header_cell = ["Lon", "Lat"] + list(out_array["varNames"])
        else:
            raise RuntimeError("in_data has neither garr_xvy nor garr_xv")
        return in_header_cell, None

    if isinstance(header, str):
        with open(header, "r") as f:
            in_header_str = f.readline().rstrip("\n")
        in_header_cell = in_header_str.split()
        return in_header_cell, in_header_str

    if isinstance(header, (list, tuple)):
        if len(header) == 2 and isinstance(header[1], (list, tuple)):
            return list(header[1]), header[0]
        return list(header), None

    raise RuntimeError("header not valid.")


def _get_lat_lon_year_cols(header_cell: Sequence[str]) -> tuple[int, int, int]:
    stripped = [h.strip() for h in header_cell]
    if "lon" in stripped:
        lon_name = "lon"
    elif "Lon" in stripped:
        lon_name = "Lon"
    elif "Lon." in stripped:
        lon_name = "Lon."
    else:
        raise RuntimeError("Could not identify longitude column!")

    if "lat" in stripped:
        lat_name = "lat"
    elif "Lat" in stripped:
        lat_name = "Lat"
    elif "Lat." in stripped:
        lat_name = "Lat."
    else:
        raise RuntimeError("Could not identify latitude column!")

    if "year" in stripped:
        year_name = "year"
    elif "Year" in stripped:
        year_name = "Year"
    else:
        year_name = ""

    i_lon = stripped.index(lon_name)
    i_lat = stripped.index(lat_name)
    i_year = stripped.index(year_name) if year_name else -1
    return i_lat, i_lon, i_year


def _import_to_table(in_file: str, verbose: bool, dispPrefix: str, force_as_gridlist: bool) -> pd.DataFrame:
    in_header = _read_header(in_file)
    if verbose:
        print(f"{dispPrefix}      Reading data...")

    with open(in_file, "r") as f:
        firstline = f.tell()
        numcols = len(f.readline().strip().split())
        if len(in_header) == 2:
            f.seek(firstline)
        data = np.genfromtxt(f, dtype=float)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    # Remove trailing all-NaN columns
    while data.shape[1] > 0 and np.all(np.isnan(data[:, -1])):
        data = data[:, :-1]
        numcols -= 1

    if data.shape[1] < 2:
        raise RuntimeError("Input must be 2-d matrix.")
    if force_as_gridlist:
        data = data[:, :2]
        in_header = in_header[:2]
        col_names = ["Lon", "Lat"]
    elif data.shape[1] == 2:
        warnings.warn("Assuming columns Lon, Lat.")
        col_names = ["Lon", "Lat"]
    else:
        col_names = in_header

    if len(in_header) != data.shape[1]:
        raise RuntimeError("Length of column name list does not match horizontal dimension of array.")

    var_names = [c.replace('"', "").replace("-", "_") for c in col_names]
    out_table = pd.DataFrame(data, columns=var_names)

    # Special case for doubled 2005 rows
    base = Path(in_file).name
    if base == "N_fert_rcp85_6f.out" and "Year" in out_table.columns:
        bad = out_table.index[out_table["Year"] == 2005].to_list()
        if bad:
            out_table = out_table.drop(bad[1::2]).reset_index(drop=True)

    return _standardize_colnames(out_table)


def _read_header(in_file: str) -> list[str]:
    with open(in_file, "r") as f:
        line = f.readline()
    header = line.split()
    if header and header[0] == "":
        header = header[1:]
    if header and header[-1] == "":
        header = header[:-1]
    for h in header:
        if h == "":
            raise RuntimeError("Empty variable name!")
    return header


def _standardize_colnames(df: pd.DataFrame) -> pd.DataFrame:
    renames = {}
    for col in df.columns:
        if col == "lon":
            renames[col] = "Lon"
        elif col == "lat":
            renames[col] = "Lat"
        elif col == "year":
            renames[col] = "Year"
    if renames:
        df = df.rename(columns=renames)
    return df


def _gunzip_if_needed(in_file: str, verbose: bool, dispPrefix: str) -> bool:
    did_unzip = False
    if not os.path.exists(in_file):
        gz_path = f"{in_file}.gz"
        gz_path = _resolve_symlink(gz_path)
        if os.path.exists(gz_path):
            if verbose:
                print(f"{dispPrefix}   Unzipping in_file...")
            with gzip.open(gz_path, "rb") as f_in, open(in_file, "wb") as f_out:
                f_out.write(f_in.read())
            did_unzip = True
        else:
            raise FileNotFoundError(f"{in_file}.mat, {in_file}, and {in_file}.gz not found.")
    return did_unzip


def _resolve_symlink(path: str) -> str:
    if os.path.islink(path):
        return os.path.realpath(path)
    return path


def _expand_wildcard(path: str) -> str:
    if "*" not in path:
        return path
    matches = sorted(Path().glob(path))
    if len(matches) == 0:
        raise FileNotFoundError(f"No match found for {path}")
    if len(matches) > 1:
        raise RuntimeError(f"More than one match found for {path}")
    return str(matches[0].resolve())


def _strip_known_suffixes(path: str) -> str:
    if path.endswith(".gz"):
        path = path[:-3]
    if path.endswith(".mat"):
        path = path[:-4]
    if path.endswith(".maps"):
        path = path[:-5]
    return path


def _get_indices(
    in_table: pd.DataFrame,
    xres: float,
    yres: float,
    list_to_map_in: Sequence[int] | None,
    lat_orient: str,
    lon_orient: str,
    drop_northpole: bool,
    drop_southpole: bool,
    lons_centered_on_180: bool,
    verboseIfNoMat: bool,
    verbose: bool,
    in_prec: int,
) -> tuple[list[int], bool, list[str], np.ndarray, float, float, np.ndarray, list[float]]:
    in_lons = in_table["Lon"].to_numpy()
    in_lats = in_table["Lat"].to_numpy()
    in_ncells = len(in_lats)

    col_names, var_names, multi_yrs, year_list = _get_names(in_table)
    if multi_yrs:
        in_years = in_table["Year"].to_numpy()
        min_year = np.min(in_years)
        in_lons = in_lons[in_years == min_year]
        in_lats = in_lats[in_years == min_year]
        in_ncells = len(in_lats)

    xres, yres, lat_extent = _process_resolution(
        xres, yres, in_lons, in_lats, lat_orient, drop_northpole, drop_southpole, verboseIfNoMat, verbose
    )

    lons_map, lats_map = _set_up_maps(
        xres,
        yres,
        in_lons,
        in_lats,
        lat_orient,
        lon_orient,
        lat_extent,
        lons_centered_on_180,
        verboseIfNoMat,
        verbose,
    )

    if list_to_map_in is None or len(list_to_map_in) == 0:
        list_to_map, found = _get_map_indices(in_lons, in_lats, lons_map, lats_map, verboseIfNoMat, verbose, in_prec)
    else:
        if len(in_lons) != len(list_to_map_in):
            warnings.warn("length(in_lons) ~= length(list_to_map_in)! Ignoring list_to_map_in.")
            list_to_map, found = _get_map_indices(
                in_lons, in_lats, lons_map, lats_map, verboseIfNoMat, verbose, in_prec
            )
        else:
            list_to_map = np.array(list_to_map_in, dtype=int)
            found = np.ones_like(list_to_map, dtype=bool)

    return year_list, multi_yrs, var_names, list_to_map, xres, yres, found, lat_extent


def _get_names(in_table: pd.DataFrame) -> tuple[list[str], list[str], bool, list[int]]:
    col_names = list(in_table.columns)
    var_names = []
    multi_yrs = False
    year_list: list[int] = []
    for name in col_names:
        if name not in {"Lat", "Lon", "Year"}:
            var_names.append(name)
        elif name == "Year":
            years = np.unique(in_table["Year"].to_numpy())
            if len(years) > 1:
                multi_yrs = True
                year_list = years.tolist()
    return col_names, var_names, multi_yrs, year_list


def _process_resolution(
    xres_in: float,
    yres_in: float,
    in_lons: np.ndarray,
    in_lats: np.ndarray,
    lat_orient: str,
    drop_northpole: bool,
    drop_southpole: bool,
    verboseIfNoMat: bool,
    verbose: bool,
) -> tuple[float, float, list[float]]:
    unique_lats = np.unique(in_lats)
    if (drop_northpole or drop_southpole) and lat_orient != "center":
        raise RuntimeError(f"You asked to drop north and/or south pole, but lat_orient is {lat_orient}")
    if drop_northpole and np.any(unique_lats == 90):
        raise RuntimeError("drop_northpole requested but 90° is in in_lats. Code to deal with that.")
    if drop_southpole and np.any(unique_lats == -90):
        raise RuntimeError("drop_southpole requested but -90° is in in_lats. Code to deal with that.")

    if xres_in > 0 and yres_in > 0:
        xres_out = xres_in
        yres_out = yres_in
    elif xres_in > 0 and not (yres_in > 0):
        xres_out = xres_in
        yres_out = xres_in
    elif not (xres_in > 0) and yres_in > 0:
        xres_out = yres_in
        yres_out = yres_in
    else:
        if verboseIfNoMat or verbose:
            print("      Determining X and Y resolution...")
        if not (yres_in > 0):
            diffs = np.abs(np.diff(unique_lats))
            yres_out = np.min(diffs)
        if not (xres_in > 0):
            unique_lons = np.unique(in_lons)
            diffs = np.abs(np.diff(unique_lons))
            xres_out = np.min(diffs)
        if verboseIfNoMat or verbose:
            print(f"      Assuming X res. = {xres_out}, Y res. = {yres_out}")

    lat_extent = _get_lat_extent(lat_orient, drop_northpole, drop_southpole, yres_out)
    return xres_out, yres_out, lat_extent


def _get_lat_extent(lat_orient: str, drop_northpole: bool, drop_southpole: bool, yres: float) -> list[float]:
    if (drop_northpole or drop_southpole) and lat_orient != "center":
        raise RuntimeError(f"You asked to drop north and/or south pole, but lat_orient is {lat_orient}")
    lat_extent = [-90.0, 90.0]
    if lat_orient == "center":
        if drop_northpole:
            lat_extent[1] = lat_extent[1] - yres / 2
        if drop_southpole:
            lat_extent[0] = lat_extent[0] + yres / 2
    elif lat_orient not in {"upper", "lower", ""} and (drop_southpole or drop_northpole):
        raise RuntimeError(f"lpjgu_process_resolution() doesn't know how to handle lat_orient {lat_orient}")
    return lat_extent


def _set_up_maps(
    xres: float,
    yres: float,
    in_lons: np.ndarray,
    in_lats: np.ndarray,
    lat_orient: str,
    lon_orient: str,
    lat_extent: Sequence[float],
    lons_centered_on_180: bool,
    verboseIfNoMat: bool,
    verbose: bool,
) -> tuple[np.ndarray, np.ndarray]:
    if lat_orient == "":
        if np.any(in_lats == lat_extent[0]):
            lat_orient = "lower"
        elif np.any(in_lats == lat_extent[1]):
            lat_orient = "upper"
        elif np.any((in_lats - yres / 2 == lat_extent[0]) | (in_lats + yres / 2 == lat_extent[1])):
            lat_orient = "center"
        else:
            if np.any(np.mod(in_lats, yres) == 0):
                lat_orient = "lower"
            else:
                lat_orient = "center"
        if verboseIfNoMat or verbose:
            print(f"      Assuming lat_orient = {lat_orient}.")

    if lon_orient == "":
        if lons_centered_on_180:
            if np.any(in_lons + xres / 2 == 180):
                lon_orient = "left"
            elif np.any(in_lons == 180):
                lon_orient = "center"
            elif np.any(in_lons - xres / 2 == -180):
                lon_orient = "right"
            else:
                raise RuntimeError("Figure this out for lons_centered_on_180.")
        else:
            if np.any(in_lons == -180):
                lon_orient = "left"
            elif np.any(in_lons == 180):
                lon_orient = "right"
            elif np.any((in_lons - xres / 2 == -180) | (in_lons + xres / 2 == 180)):
                lon_orient = "center"
            else:
                if np.any(np.mod(in_lons, xres) == 0):
                    lon_orient = "left"
                else:
                    lon_orient = "center"
        if verboseIfNoMat or verbose:
            print(f"      Assuming lon_orient = {lon_orient}.")

    if lat_orient == "lower":
        lat_min = lat_extent[0]
        lat_max = lat_extent[1] - yres
    elif lat_orient == "upper":
        lat_min = lat_extent[0] + yres
        lat_max = lat_extent[1]
    else:
        lat_min = lat_extent[0] + yres / 2
        lat_max = lat_extent[1] - yres / 2

    if lons_centered_on_180:
        if lon_orient in {"left", "right"}:
            raise RuntimeError("Figure this out for lons_centered_on_180 left/right.")
        if np.any(in_lons == 180):
            lon_min = -180 + xres
            lon_max = 180
        else:
            lon_min = -180
            lon_max = 180 - xres
    else:
        if lon_orient == "left":
            lon_min = -180
            lon_max = 180 - xres
        elif lon_orient == "right":
            lon_min = -180 + xres
            lon_max = 180
        else:
            lon_min = -180 + xres / 2
            lon_max = 180 - xres / 2

    lons = np.arange(lon_min, lon_max + 1e-9, xres)
    lats = np.arange(lat_min, lat_max + 1e-9, yres)
    out_lons_map = np.tile(lons, (len(lats), 1))
    out_lats_map = np.tile(lats.reshape(-1, 1), (1, len(lons)))
    return out_lons_map, out_lats_map


def _get_map_indices(
    in_lons: np.ndarray,
    in_lats: np.ndarray,
    in_lons_map: np.ndarray,
    in_lats_map: np.ndarray,
    verboseIfNoMat: bool,
    verbose: bool,
    in_prec: int,
) -> tuple[np.ndarray, np.ndarray]:
    if verboseIfNoMat or verbose:
        print("      Getting indices to convert list to map...")
    n_cells = len(in_lons)
    lons_vec = in_lons_map[0, :]
    lats_vec = in_lats_map[:, 0]

    in_lons = np.round(in_lons, in_prec)
    in_lats = np.round(in_lats, in_prec)
    lons_vec = np.round(lons_vec, in_prec)
    lats_vec = np.round(lats_vec, in_prec)

    lon_index_map = {val: idx for idx, val in enumerate(lons_vec)}
    lat_index_map = {val: idx for idx, val in enumerate(lats_vec)}
    lons_inds = np.array([lon_index_map.get(x, -1) for x in in_lons])
    lats_inds = np.array([lat_index_map.get(x, -1) for x in in_lats])
    found = (lons_inds >= 0) & (lats_inds >= 0)

    if not np.any(found):
        raise RuntimeError(
            f"No cells found with both lon and lat! {np.sum(lons_inds>=0)} found lon, {np.sum(lats_inds>=0)} found lat"
        )
    if np.any(~found):
        warnings.warn(f"{np.sum(~found)} cells being ignored.")

    lons_inds = lons_inds[found]
    lats_inds = lats_inds[found]
    list_to_map = np.ravel_multi_index((lats_inds, lons_inds), in_lons_map.shape, order="F") + 1

    if np.any(np.isnan(list_to_map)):
        raise RuntimeError("Somehow list_to_map contains NaN.")
    if len(list_to_map) != n_cells - np.sum(~found):
        raise RuntimeError("length(list_to_map) ~= Ncells-length(find(~found))")
    if len(list_to_map) != len(np.unique(list_to_map)):
        raise RuntimeError("list_to_map has duplicates")

    return list_to_map, found


def _get_nlonlat_in_tolerance(n_deg: float, res: float, quiet: bool) -> int:
    n = n_deg / res
    if round(n) != n:
        tol = 1e-3
        if abs(round(n) - n) > tol:
            raise RuntimeError(
                f"Nlon/lat not integer within {tol}: off from {round(n)} by {abs(round(n)-n)} (Ndeg {n_deg}, res {res})"
            )
        if not quiet:
            warnings.warn(f"Nlon/lat not integer (off by {abs(round(n)-n)}) but within tolerance of {tol}")
    return int(round(n))


def _make_maps(
    xres: float,
    yres: float,
    multi_yrs: bool,
    yearList: Sequence[int],
    varNames: Sequence[str],
    in_table: pd.DataFrame,
    list_to_map: Sequence[int],
    is_gridlist: bool,
    found: np.ndarray,
    lat_extent: Sequence[float],
    dataType: str,
    verboseIfNoMat: bool,
    verbose: bool,
    quiet: bool,
):
    nlon = _get_nlonlat_in_tolerance(360, xres, quiet)
    nlat = _get_nlonlat_in_tolerance(lat_extent[1] - lat_extent[0], yres, quiet)
    if is_gridlist:
        if verboseIfNoMat or verbose:
            print("      Making mask...")
        out_maps = np.zeros((nlat, nlon), dtype=bool)
        rr, cc = np.unravel_index(list_to_map - 1, out_maps.shape, order="F")
        out_maps[rr, cc] = True
        return out_maps

    if verboseIfNoMat or verbose:
        print("      Making maps...")
    nvars = len(varNames)
    dtype = np.float64 if dataType == "double" else np.float32
    if multi_yrs:
        nyears = len(yearList)
        out_maps = np.full((nlat, nlon, nvars, nyears), np.nan, dtype=dtype)
    else:
        out_maps = np.full((nlat, nlon, nvars), np.nan, dtype=dtype)

    for v, this_var in enumerate(varNames):
        if verboseIfNoMat or verbose:
            print(f"         {this_var} ({v+1} of {nvars})...")
        if multi_yrs:
            for y, this_year in enumerate(yearList):
                tmp = np.full((nlat, nlon), np.nan, dtype=dtype)
                vals = in_table.loc[in_table["Year"] == this_year, this_var].to_numpy()
                rr, cc = np.unravel_index(list_to_map - 1, tmp.shape, order="F")
                tmp[rr, cc] = vals
                out_maps[:, :, v, y] = tmp
        else:
            tmp = np.full((nlat, nlon), np.nan, dtype=dtype)
            vals = in_table.loc[found, this_var].to_numpy()
            rr, cc = np.unravel_index(list_to_map - 1, tmp.shape, order="F")
            tmp[rr, cc] = vals
            out_maps[:, :, v] = tmp

    return out_maps


def _maps_to_garr(in_struct: dict) -> dict:
    if "maps_YXvy" in in_struct:
        maps = in_struct["maps_YXvy"]
        nvars = maps.shape[2]
        mask = ~np.isnan(np.mean(maps, axis=(2, 3)))
        list2map = np.flatnonzero(mask.ravel(order="F")) + 1
        list2map_all = _get_list2map_all(list2map, maps.shape[:2], maps.shape[2:])
        tmp_xvy = np.full((len(list2map),) + maps.shape[2:], np.nan)
        tmp_xvy.ravel(order="F")[:] = maps.ravel(order="F")[list2map_all - 1]
        out = {
            "list2map": list2map,
            "varNames": in_struct.get("varNames", []),
            "map_size": maps.shape[:2],
            "yearList": in_struct.get("yearList", []),
            "garr_xvy": tmp_xvy,
        }
    elif "maps_YXv" in in_struct:
        maps = in_struct["maps_YXv"]
        mask = ~np.isnan(np.mean(maps, axis=2))
        list2map = np.flatnonzero(mask.ravel(order="F")) + 1
        list2map_all = _get_list2map_all(list2map, maps.shape[:2], maps.shape[2:])
        tmp_xv = np.full((len(list2map), maps.shape[2]), np.nan)
        tmp_xv.ravel(order="F")[:] = maps.ravel(order="F")[list2map_all - 1]
        out = {
            "list2map": list2map,
            "varNames": in_struct.get("varNames", []),
            "map_size": maps.shape[:2],
            "garr_xv": tmp_xv,
        }
    else:
        raise RuntimeError("maps_YXv(y) not found")
    if "lat_extent" in in_struct:
        out["lat_extent"] = in_struct["lat_extent"]
    if "lat_orient" in in_struct:
        out["lat_orient"] = in_struct["lat_orient"]
    if "lonlats" in in_struct:
        out["lonlats"] = in_struct["lonlats"]
    return out


def _get_list2map_all(list2map: np.ndarray, map_size: Sequence[int], rep_sizes: Sequence[int]) -> np.ndarray:
    n_reps = int(np.prod(rep_sizes))
    n_cells = len(list2map)
    list2map_all = np.empty(n_cells * n_reps, dtype=float)
    for r in range(n_reps):
        i1 = r * n_cells
        iN = (r + 1) * n_cells
        list2map_all[i1:iN] = r * (map_size[0] * map_size[1]) + list2map
    if np.any(np.isnan(list2map_all)):
        raise RuntimeError("NaN(s) in list2map_all")
    if np.any(list2map_all <= 0):
        raise RuntimeError("Non-positive value(s) in list2map_all")
    if len(np.unique(list2map_all)) != len(list2map_all):
        raise RuntimeError("Non-unique(s) in list2map_all")
    return list2map_all.astype(int)


def _mat_struct_to_dict(mat_obj) -> dict:
    # Minimal conversion for scipy.io.loadmat outputs.
    if isinstance(mat_obj, np.ndarray) and mat_obj.dtype.names:
        out = {}
        for name in mat_obj.dtype.names:
            out[name] = mat_obj[name][0, 0]
        return out
    if isinstance(mat_obj, dict):
        return mat_obj
    return {"out_struct": mat_obj}


def _is_hdf5_mat(path: str) -> bool:
    try:
        with h5py.File(path, "r"):
            return True
    except Exception:
        return False


def _read_maps_mat_h5(path: str) -> dict:
    out = {}
    with h5py.File(path, "r") as f:
        if "out_struct" not in f:
            return {}
        g = f["out_struct"]
        if "maps_YXv" in g:
            maps = np.array(g["maps_YXv"])
            # MATLAB stores as YXv; h5py often returns vX Y
            if "varNames" in g:
                vnames = _read_cellstr_h5(f, g["varNames"])
                if maps.ndim == 3 and maps.shape[0] == len(vnames):
                    maps = np.transpose(maps, (2, 1, 0))
                    out["varNames"] = vnames
            out["maps_YXv"] = maps
        if "maps_YXvy" in g:
            maps = np.array(g["maps_YXvy"])
            if "varNames" in g:
                vnames = _read_cellstr_h5(f, g["varNames"])
                if maps.ndim == 4 and maps.shape[0] == len(vnames):
                    maps = np.transpose(maps, (3, 2, 1, 0))
                    out["varNames"] = vnames
            out["maps_YXvy"] = maps
        if "varNames" in g and "varNames" not in out:
            out["varNames"] = _read_cellstr_h5(f, g["varNames"])
        if "yearList" in g:
            out["yearList"] = np.array(g["yearList"]).ravel()
        if "list_to_map" in g:
            out["list_to_map"] = np.array(g["list_to_map"]).ravel()
        if "lonlats" in g:
            out["lonlats"] = np.array(g["lonlats"])
        if "mask_YX" in g:
            out["mask_YX"] = np.array(g["mask_YX"])
        if "lat_extent" in g:
            out["lat_extent"] = np.array(g["lat_extent"]).ravel()
        if "lat_orient" in g:
            out["lat_orient"] = _read_char_h5(f, g["lat_orient"])
    return out


def _read_cellstr_h5(root: h5py.File, ds: h5py.Dataset) -> list[str]:
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


def _normalize_garr_struct(out_struct: dict) -> dict:
    # Handle garr_yxv -> garr_xvy conversion and remove single-year dimension.
    if "garr_yxv" in out_struct:
        tmp = np.transpose(out_struct["garr_yxv"], (1, 2, 0))
        out_struct.pop("garr_yxv")
        out_struct["garr_xvy"] = tmp
    if "garr_xvy" in out_struct and out_struct["garr_xvy"].shape[2] == 1:
        tmp = out_struct["garr_xvy"]
        out_struct.pop("garr_xvy")
        out_struct["garr_xv"] = tmp
    return out_struct


def _table_to_garr_with_target(table_in: pd.DataFrame, target) -> dict:
    if isinstance(target, dict):
        lonlats_target = target.get("lonlats")
        list2map_target = target.get("list2map")
    elif isinstance(target, (list, tuple)) and len(target) == 2:
        lonlats_target, list2map_target = target
    else:
        raise RuntimeError("Input target is of invalid class")
    if lonlats_target is None or list2map_target is None:
        raise RuntimeError("target must provide lonlats and list2map")

    lonlats_target = np.asarray(lonlats_target)
    list2map_target = np.asarray(list2map_target, dtype=int).ravel()

    lonlats_in = table_in[["Lon", "Lat"]].drop_duplicates().to_numpy()
    n_cells = lonlats_in.shape[0]

    year_list = np.unique(table_in["Year"].to_numpy())
    n_years = len(year_list)

    var_names = [c for c in table_in.columns if c not in {"Lon", "Lat", "Year"}]
    n_vars = len(var_names)

    if len(table_in) != n_cells * n_years:
        raise RuntimeError("length(table_in.Lon) ~= Ncells*Nyears")

    vals = table_in[var_names].to_numpy()
    garr_yxv = np.empty((n_years, n_cells, n_vars))
    for v in range(n_vars):
        garr_yxv[:, :, v] = vals[:, v].reshape((n_years, n_cells), order="F")
    garr_xvy = np.transpose(garr_yxv, (1, 2, 0))

    # Reorder to match target lonlats
    if not np.array_equal(lonlats_target, lonlats_in):
        index_map = {tuple(row): idx for idx, row in enumerate(lonlats_in)}
        if len(index_map) != len(lonlats_in):
            raise RuntimeError("Duplicate lon/lat rows in input table")
        idxs = []
        for row in lonlats_target:
            key = (row[0], row[1])
            if key not in index_map:
                raise RuntimeError(
                    f"Gridlist mismatch: lon/lat {row[0]:0.2f} {row[1]:0.2f} missing from input table"
                )
            idxs.append(index_map[key])
        garr_xvy = garr_xvy[np.array(idxs, dtype=int), :, :]

    return {
        "list2map": list2map_target,
        "lonlats": lonlats_target,
        "varNames": var_names,
        "yearList": year_list,
        "garr_xvy": garr_xvy,
    }


def _load_table_from_mat(in_matfile: str) -> pd.DataFrame:
    loaded = loadmat(in_matfile)
    if "out_table" not in loaded:
        raise RuntimeError("MAT file missing out_table")
    out = loaded["out_table"]
    if isinstance(out, np.ndarray) and out.dtype.names:
        data = {name: out[name].squeeze() for name in out.dtype.names}
        return pd.DataFrame(data)
    if isinstance(out, np.ndarray):
        return pd.DataFrame(out)
    raise RuntimeError("Unsupported MAT table format")
