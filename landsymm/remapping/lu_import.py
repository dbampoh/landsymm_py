"""Land-use import utilities (stub).

Pseudocode:
- import_lu_hildaplus:
  1) Read HILDA+ netfrac table to geoarray (target gridlist).
  2) Subset to yearList_out.
  3) Map input classes to NATURAL/CROPLAND/PASTURE/BARREN.
  4) Add water fraction to BARREN from staticData_quarterdeg.nc.
  5) Validate sum-to-1.

- import_lu_luh2:
  1) Read LUH2 NetCDF states.
  2) Map input classes to NATURAL/CROPLAND/PASTURE/BARREN.
  3) Aggregate to 0.5-degree using area weights.
  4) Add water fraction to BARREN.

- note_unveg:
  1) Identify fully BARREN cells.
  2) If fill_unveg < 0 -> error.
  3) If fill_unveg > 0 -> transfer area from BARREN to NATURAL.
"""

# TODO:
# - implement HILDA+ and LUH2 import with exact class mapping
# - match MATLAB aggregation and water-fraction handling
from __future__ import annotations

import os
from typing import Sequence, Tuple

import numpy as np
from netCDF4 import Dataset

from landsymm.common.aggregation import aggregate
from landsymm.common.geoarray import GeoArray
from landsymm.common.lpjg_io import read2geoarray
from landsymm.common.mapping_tools import yxz_to_xz


def import_lu_hildaplus(
    geodata_dir: str, year_list_out: Sequence[int], gridlist: GeoArray
) -> Tuple[GeoArray, np.ndarray]:
    """Import HILDA+ netfrac and map to LU classes (TODO).

    Pseudocode:
    1) Read netfrac table with target gridlist.
    2) Subset to year_list_out.
    3) Map input LU classes to NATURAL/CROPLAND/PASTURE/BARREN.
    4) Add ice/water fraction to BARREN.
    5) Validate sum-to-1.
    """
    year_list_lu_states = np.arange(1901, 2021)
    file_lu_states = (
        f"{geodata_dir}/HILDA+/data/output/"
        "hildaplus_netfrac_1901_2020.txt"
    )
    if min(year_list_out) < year_list_lu_states.min() or max(year_list_out) > year_list_lu_states.max():
        raise RuntimeError(
            "year_list_out (%d-%d) not entirely contained within yearList_lu_states (%d-%d)"
            % (min(year_list_out), max(year_list_out), year_list_lu_states.min(), year_list_lu_states.max())
        )

    print("Importing land uses...")
    print(f"  file_lu_states: {file_lu_states}")
    in_lu = read2geoarray(
        file_lu_states,
        target=gridlist,
        verboseIfNoMat=False,
        force_mat_save=True,
        force_mat_nosave=True,
    )
    # Match MATLAB text parsing precision for HILDA+ inputs
    if "garr_xvy" in in_lu:
        in_lu["garr_xvy"] = np.round(in_lu["garr_xvy"], 7)

    # Subset years
    year_list_out = list(year_list_out)
    year_list_in = np.ravel(np.asarray(in_lu["yearList"])).tolist()
    year_list_in = [int(y) for y in year_list_in]
    ib = [year_list_in.index(int(y)) for y in year_list_out]
    if [year_list_in[i] for i in ib] != year_list_out:
        raise RuntimeError("yearList mismatch")
    in_lu["garr_xvy"] = in_lu["garr_xvy"][:, :, ib]
    in_lu["yearList"] = np.array(year_list_out)

    list_lu_in = list(in_lu["varNames"])
    list_lu_out = ["NATURAL", "CROPLAND", "PASTURE", "BARREN"]

    map_lu_in2out = []
    for name in list_lu_in:
        if "CROPLAND" in name:
            map_lu_in2out.append("CROPLAND")
        elif "PASTURE" in name:
            map_lu_in2out.append("PASTURE")
        elif ("NATURAL" in name) or ("FOREST" in name):
            map_lu_in2out.append("NATURAL")
        elif ("URBAN" in name) or ("BARREN" in name):
            map_lu_in2out.append("BARREN")
        else:
            map_lu_in2out.append("")

    if any(m == "" for m in map_lu_in2out):
        raise RuntimeError("Remaining empty element(s) in map_LU_in2out!")
    if len(set(map_lu_in2out) - set(list_lu_out)) > 0:
        raise RuntimeError("Some member of map_LU_in2out is not present in list_LU_out!")
    print("Warning: Should work out specification of mapping with check for duplicates on LHS.")

    if np.isnan(in_lu["garr_xvy"]).any():
        raise RuntimeError("Handle NaNs in in_lu.garr_xvy")

    n_cells = len(gridlist["list2map"])
    n_years = len(year_list_out)
    out_lu = {k: v for k, v in in_lu.items() if k not in {"garr_xvy", "varNames"}}
    out_lu["varNames"] = list_lu_out
    out_lu["garr_xvy"] = np.zeros((n_cells, len(list_lu_out), n_years))
    for v, this_in in enumerate(list_lu_in):
        this_out = map_lu_in2out[v]
        print(f"    {this_in} ({v+1} of {len(list_lu_in)}) to {this_out}...")
        i = list_lu_out.index(this_out)
        out_lu["garr_xvy"][:, i, :] += in_lu["garr_xvy"][:, v, :]

    if _max_sum_diff(out_lu["garr_xvy"], 1.0) > 1e-6:
        raise RuntimeError("Combining input LU types into output types: Changed sum of LU fractions")

    # Import cell area and add water fraction
    file_lu_etc = _find_staticdata_path(geodata_dir)
    print(f"file_lu_etc: {file_lu_etc}")
    with Dataset(file_lu_etc, "r") as ds:
        carea_XY = np.asarray(ds.variables["carea"][:], dtype=float)
        icwtr_XY = np.asarray(ds.variables["icwtr"][:], dtype=float)

    # netCDF arrays are lat x lon; MATLAB expects lon x lat then transposes+flipud.
    carea_YX = np.flipud(carea_XY)
    carea_YX = aggregate(carea_YX, 0.25, 0.5)
    carea_hd_YX = aggregate(carea_XY, 0.25, 0.5)

    icwtr_YX = np.flipud(icwtr_XY)
    icwtr_YX[icwtr_YX == 1] = 0
    icwtr_hd_YX = aggregate(icwtr_YX * np.flipud(carea_XY), 0.25, 0.5) / np.flipud(carea_hd_YX)
    icwtr_hd_x = icwtr_hd_YX.ravel(order="F")[np.asarray(gridlist["list2map"]) - 1]
    v_bare = list_lu_out.index("BARREN")
    icwtr_hd_x1y = np.tile(icwtr_hd_x.reshape(-1, 1), (1, n_years))
    icwtr_hd_xvy = np.tile(icwtr_hd_x1y.reshape(n_cells, 1, n_years), (1, len(list_lu_out), 1))

    out_lu["garr_xvy"] = out_lu["garr_xvy"] * (1 - icwtr_hd_xvy)
    out_lu["garr_xvy"][:, v_bare, :] = out_lu["garr_xvy"][:, v_bare, :] + icwtr_hd_x1y

    if _max_sum_diff(out_lu["garr_xvy"], 1.0) >= 1e-6:
        raise RuntimeError("Error in adding ice/water fraction to BARREN: Fractions do not sum to 1")

    print("Done.")
    return out_lu, carea_YX


def import_lu_luh2(
    geodata_dir: str, year_list_out: Sequence[int], gridlist: GeoArray
) -> Tuple[GeoArray, np.ndarray]:
    """Import LUH2 netcdf and map to LU classes (TODO)."""
    year_list_out = list(year_list_out)
    year_list_lu_states = list(range(1850, 2016))
    file_lu_states = f"{geodata_dir}/LUH2/v2h/states.{year_list_lu_states[0]}-{year_list_lu_states[-1]}.nc"
    if min(year_list_out) < year_list_lu_states[0] or max(year_list_out) > year_list_lu_states[-1]:
        raise RuntimeError("year_list_out must be entirely contained within yearList_lu_states!")

    year_is = [year_list_lu_states.index(y) for y in year_list_out]
    n_years = len(year_list_out)

    with Dataset(file_lu_states, "r") as ds:
        list_lu_in = [v for v in ds.variables.keys() if v not in {"time", "lon", "lat", "lat_bounds", "lon_bounds", "secma", "secmb"}]

        list_lu_out = ["NATURAL", "CROPLAND", "PASTURE", "BARREN"]
        map_lu_in2out = []
        for name in list_lu_in:
            if any(k in name for k in ["c3ann", "c3nfx", "c3per", "c4ann", "c4per"]):
                map_lu_in2out.append("CROPLAND")
            elif any(k in name for k in ["pastr", "range"]):
                map_lu_in2out.append("PASTURE")
            elif any(k in name for k in ["primf", "primn", "secdf", "secdn", "secma", "secmb"]):
                map_lu_in2out.append("NATURAL")
            elif "urban" in name:
                map_lu_in2out.append("BARREN")
            else:
                map_lu_in2out.append("")

        if any(m == "" for m in map_lu_in2out):
            raise RuntimeError("Remaining empty element(s) in map_LU_in2out!")
        if len(set(map_lu_in2out) - set(list_lu_out)) > 0:
            raise RuntimeError("Some member of map_LU_in2out is not present in list_LU_out!")
        print("Warning: Should work out specification of mapping with check for duplicates on LHS.")

        print("Importing land uses...")
        print(f"  file_lu_states: {file_lu_states}")

        file_lu_etc = _find_staticdata_path(geodata_dir)
        print(f"file_lu_etc: {file_lu_etc}")
        with Dataset(file_lu_etc, "r") as ds_etc:
            carea_XY = np.asarray(ds_etc.variables["carea"][:], dtype=float)
            icwtr_XY = np.asarray(ds_etc.variables["icwtr"][:], dtype=float)

        carea_YX = np.flipud(carea_XY)
        carea_YX = aggregate(carea_YX, 0.25, 0.5)
        carea_XYy = np.repeat(carea_XY[:, :, None], n_years, axis=2)
        carea_hd_YX = aggregate(carea_XY, 0.25, 0.5)
        carea_hd_YXy = np.repeat(carea_hd_YX[:, :, None], n_years, axis=2)

        out_lu = {k: v for k, v in gridlist.items() if k != "mask_YX"}
        out_lu["varNames"] = list_lu_out
        out_lu["yearList"] = np.array(year_list_out)
        n_cells = len(gridlist["list2map"])

        for v, this_in in enumerate(list_lu_in):
            this_out = map_lu_in2out[v]
            print(f"    {this_in} ({v+1} of {len(list_lu_in)}) to {this_out}...")
            i = list_lu_out.index(this_out)

            lu_in_XYy = ds.variables[this_in][:, :, year_is]
            lu_in_XYy = np.nan_to_num(lu_in_XYy, nan=0.0)
            lu_out_YXy = np.flipud(
                np.transpose(
                    aggregate(lu_in_XYy * carea_XYy, 0.25, 0.5) / carea_hd_XYy, (1, 0, 2)
                )
            )
            if (lu_out_YXy - 1 > 1e-6).any():
                raise RuntimeError("Some element(s) of lu_out_YXy > 1!")

            lu_out_xy = yxz_to_xz(lu_out_YXy, (n_cells, lu_out_YXy.shape[2]), gridlist["list2map"])
            if v == 0:
                out_lu["garr_xvy"] = np.zeros((n_cells, len(list_lu_out), lu_out_xy.shape[1]))
            out_lu["garr_xvy"][:, i, :] += np.transpose(lu_out_xy.reshape(n_cells, 1, -1), (0, 1, 2))
            if (out_lu["garr_xvy"] - 1 > 1e-6).any():
                raise RuntimeError("Some element(s) of out_lu.garr_xvy > 1!")

        icwtr_YX = np.flipud(icwtr_XY)
        icwtr_YX[icwtr_YX == 1] = 0
        icwtr_hd_YX = aggregate(icwtr_YX * np.flipud(carea_XY), 0.25, 0.5) / np.flipud(carea_hd_YX)
        icwtr_hd_x = icwtr_hd_YX.ravel(order="F")[np.asarray(gridlist["list2map"]) - 1]
        v_bare = list_lu_out.index("BARREN")
        out_lu["garr_xvy"][:, v_bare, :] = out_lu["garr_xvy"][:, v_bare, :] + np.tile(
            icwtr_hd_x.reshape(-1, 1), (1, n_years)
        )

    print("Done.")
    return out_lu, carea_YX


def note_unveg(out_lu: GeoArray, fill_unveg: float) -> GeoArray:
    """Handle unvegetated cells per fill_unveg (TODO)."""
    list_names = list(out_lu["varNames"])
    if "BARREN" not in list_names:
        raise RuntimeError("Expected 1 match of 'BARREN' in out_lu.varNames")
    i_bare = list_names.index("BARREN")
    sum_x1y = np.sum(out_lu["garr_xvy"][:, i_bare : i_bare + 1, :], axis=1)
    bad_x = np.max(sum_x1y, axis=1) == 1
    n_bad = int(np.sum(bad_x))
    if n_bad:
        x = int(np.where(bad_x)[0][0])
        lonlat = out_lu["lonlats"][x, :]
        msg = f"{n_bad} cells not vegetated; e.g. lon {lonlat[0]:0.2f} lat {lonlat[1]:0.2f}"

        if fill_unveg < 0:
            raise RuntimeError(msg)
        if fill_unveg > 1:
            raise RuntimeError(f"Maximum fill_unveg is 1; got {fill_unveg}")

        print(msg)
        if fill_unveg > 0:
            if "NATURAL" not in list_names:
                raise RuntimeError("Expected 1 match of 'NATURAL' in out_lu.varNames")
            i_nat = list_names.index("NATURAL")
            ntrl_x1y = out_lu["garr_xvy"][:, i_nat : i_nat + 1, :]
            bare_x1y = out_lu["garr_xvy"][:, i_bare : i_bare + 1, :]
            ntrl_x1y[bare_x1y == 1] = fill_unveg
            bare_x1y[bare_x1y == 1] = 1 - fill_unveg
            out_lu["garr_xvy"][:, i_nat : i_nat + 1, :] = ntrl_x1y
            out_lu["garr_xvy"][:, i_bare : i_bare + 1, :] = bare_x1y

            # Check that this worked
            note_unveg(out_lu, -1)

    return out_lu


def _max_sum_diff(garr_xvy: np.ndarray, target: float) -> float:
    return float(np.max(np.abs(target - np.sum(garr_xvy, axis=1))))


def _find_staticdata_path(geodata_dir: str) -> str:
    candidates = [
        f"{geodata_dir}/staticData_quarterdeg.nc",
        f"{geodata_dir}/geodata/staticData_quarterdeg.nc",
        f"{geodata_dir}/Geodata/staticData_quarterdeg.nc",
    ]
    for cand in candidates:
        if os.path.exists(cand):
            return cand
    # fallback: walk up one level
    parent = os.path.abspath(os.path.join(geodata_dir, os.pardir))
    for name in ["geodata", "Geodata"]:
        cand = os.path.join(parent, name, "staticData_quarterdeg.nc")
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError("staticData_quarterdeg.nc not found near geodata_dir")
