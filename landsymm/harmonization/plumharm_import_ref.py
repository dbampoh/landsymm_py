"""Import reference/baseline data for PLUM harmonization (MATLAB parity)."""
from __future__ import annotations

import os
from typing import Any, Dict, Sequence, Tuple

import numpy as np
from netCDF4 import Dataset

from landsymm.common.aggregation import aggregate
from landsymm.common.interpolation import inpaint_nans
from landsymm.common.lpjg_io import read_table, read_table_then2map

from .plumharm_options import PlumHarmConfig


def import_ref_data(
    cfg: PlumHarmConfig,
    *,
    do_harm: bool = True,
    year_list_baseline_lu_to_plot: Sequence[int] | None = None,
) -> Dict[str, Any]:
    """Load baseline LU/cropfracs/nfert and build mappings.

    Mirrors PLUMharm_importRefData.m.
    """
    print("Importing reference data...")

    landarea_file = os.path.join(cfg.geodata_dir, "staticData_quarterdeg.nc")
    plum_file_res_terr = os.path.join(cfg.geodata_dir, "maxcropfrac2.txt")
    plum_file_res_prot = os.path.join(cfg.geodata_dir, "protected_areas_with_points.txt")

    combine_crops = cfg.combine_crops

    # Lower-left lat/lon maps (PLUM style)
    print("    Make lower-left lat/lon map (for compat. with PLUM style)")
    lons_map_2deg = np.tile(np.arange(-180, 180, 2), (90, 1))
    lats_map_2deg = np.tile(np.arange(-90, 90, 2).reshape(-1, 1), (1, 180))
    lons_map = np.tile(np.arange(-180, 180, 0.5), (360, 1))
    lats_map = np.tile(np.arange(-90, 90, 0.5).reshape(-1, 1), (1, 720))

    # Terrain protected fraction
    print("    Read fraction of VEGETATED protected by terrain")
    res_frac_terr = read_table_then2map(
        plum_file_res_terr, verboseIfNoMat=False, force_mat_nosave=True, force_mat_save=False
    )
    res_frac_terr_yx = 1 - res_frac_terr["maps_YXv"][:, :, 0]
    res_frac_terr_yx[np.isnan(res_frac_terr_yx)] = 0

    # Protected areas
    print("    Read fraction of VEGETATED protected by protected areas")
    res_frac_prot_yx = np.loadtxt(plum_file_res_prot, skiprows=6)
    res_frac_prot_yx = np.flipud(res_frac_prot_yx)
    res_frac_prot_yx[res_frac_prot_yx == -9999] = 0

    # Land area
    print("    Get land area")
    with Dataset(landarea_file, "r") as ds:
        gcel_area_yx_qd = 1e6 * np.asarray(ds.variables["carea"][:], dtype=float)
        land_frac_yx_qd = 1 - np.flipud(np.asarray(ds.variables["icwtr"][:], dtype=float))
    land_area_yx_qd = gcel_area_yx_qd * land_frac_yx_qd

    # Convert 0.25-degree to 0.5-degree
    tmp = gcel_area_yx_qd[:, 0:1440:2] + gcel_area_yx_qd[:, 1:1440:2]
    gcel_area_yx = tmp[0:720:2, :] + tmp[1:720:2, :]
    tmp = land_area_yx_qd[:, 0:1440:2] + land_area_yx_qd[:, 1:1440:2]
    land_area_yx = tmp[0:720:2, :] + tmp[1:720:2, :]

    # Import baseline LU
    print("    Import baseline LU base_year")
    base = read_table_then2map(
        cfg.remap_lu_file, force_mat_nosave=True, force_mat_save=False
    )
    if np.any(base["maps_YXvy"][:, :, np.isin(base["varNames"], ["URBAN", "PEATLAND"]), :] > 0):
        raise RuntimeError(
            "This code is not designed to handle baseline LU inputs with any URBAN or PEATLAND area!"
        )

    if do_harm:
        tmp = base["maps_YXvy"][:, :, ~np.isin(base["varNames"], ["URBAN", "PEATLAND"]), :]
        tmp = tmp[:, :, :, base["yearList"] == cfg.base_year]
        base = {k: v for k, v in base.items() if k != "maps_YXvy"}
        base["varNames"] = [v for v in base["varNames"] if v not in {"URBAN", "PEATLAND"}]
        base["maps_YXv"] = tmp.squeeze(axis=3)
        bad_base_yx = np.sum(base["maps_YXv"], axis=2) == 0
        bad_base_yx |= np.isnan(np.sum(base["maps_YXv"], axis=2))
    else:
        if year_list_baseline_lu_to_plot is not None:
            year_list_baseline_lu_to_plot = list(year_list_baseline_lu_to_plot)
            year_is = [list(base["yearList"]).index(y) for y in year_list_baseline_lu_to_plot]
            base["maps_YXvy"] = base["maps_YXvy"][
                :, :, ~np.isin(base["varNames"], ["URBAN", "PEATLAND"]), :
            ][:, :, :, year_is]
            base["yearList"] = np.array(year_list_baseline_lu_to_plot)
        else:
            year_list_baseline_lu_to_plot = list(base["yearList"])
        base["varNames"] = [v for v in base["varNames"] if v not in {"URBAN", "PEATLAND"}]
        bad_base_yx = np.sum(base["maps_YXvy"], axis=2) == 0
        bad_base_yx |= np.isnan(np.sum(base["maps_YXvy"], axis=2))

    # Rearrange LU names
    print("    Rearrange LU names")
    lu_names = ["CROPLAND", "PASTURE", "NATURAL", "BARREN"]
    if base["varNames"] != lu_names:
        if sorted(base["varNames"]) != sorted(lu_names):
            raise RuntimeError("LUnames and base.varNames are incompatible?")
        order = [base["varNames"].index(v) for v in lu_names]
        if do_harm:
            base["maps_YXv"] = base["maps_YXv"][:, :, order]
        else:
            base["maps_YXvy"] = base["maps_YXvy"][:, :, order, :]
        base["varNames"] = lu_names
        if base["varNames"] != lu_names:
            raise RuntimeError("Error in rearranging LU names on baseline import!")

    # Harmonize masks using first PLUM dir
    print("    Harmonize masks")
    dir1 = cfg.plum_dirs[0]
    if not os.path.isdir(dir1):
        raise RuntimeError(
            f"plumDirs[0] {dir1} not found. Try changing working directory to its parent."
        )
    dir1_base = os.path.join(dir1, str(cfg.base_year))
    if not os.path.isdir(dir1_base):
        raise RuntimeError(f"Directory for {cfg.base_year} not found in {dir1_base}")
    file_in = os.path.join(dir1_base, "LandCoverFract.txt")
    S = read_table_then2map(
        file_in, verboseIfNoMat=False, force_mat_nosave=True, force_mat_save=False
    )
    mask_yx = np.isnan(S["maps_YXv"][:, :, 0])
    mask_yx |= (
        np.sum(S["maps_YXv"][:, :, np.isin(S["varNames"], ["CROPLAND", "PASTURE", "NATURAL"])], axis=2)
        == 0
    )
    mask_yx |= land_area_yx == 0
    mask_yx |= bad_base_yx
    land_area_yx[mask_yx] = 0
    if "maps_YXv" in base:
        base["maps_YXv"][np.repeat(mask_yx[:, :, None], len(lu_names), axis=2)] = np.nan
    if "maps_YXvy" in base:
        base["maps_YXvy"][
            np.repeat(mask_yx[:, :, None, None], len(lu_names), axis=2)
        ] = np.nan

    # Land area repeated per LU
    print("    Get repmat 0.5-degree land area")
    n_lu = len(lu_names)
    land_area_yxv = np.repeat(land_area_yx[:, :, None], n_lu, axis=2)

    # Avoid sum(LU) > 1 to high precision
    print("    Avoid, to high precision, sum(LU) > 1")
    eps = np.finfo(float).eps
    if do_harm:
        tmp_yx = np.sum(base["maps_YXv"], axis=2)
        j = 0
        while np.any(tmp_yx - 1 > 3 * eps):
            j += 1
            if j > 50:
                raise RuntimeError('Possible infinite loop in "Avoid, to high precision, sum(LU) > 1"')
            base["maps_YXv"] = base["maps_YXv"] / np.repeat(tmp_yx[:, :, None], n_lu, axis=2)
            tmp_yx = np.sum(base["maps_YXv"], axis=2)
    else:
        tmp_yx1y = np.sum(base["maps_YXvy"], axis=2)
        j = 0
        while np.any(tmp_yx1y - 1 > 3 * eps):
            j += 1
            if j > 50:
                raise RuntimeError('Possible infinite loop in "Avoid, to high precision, sum(LU) > 1"')
            base["maps_YXvy"] = base["maps_YXvy"] / np.repeat(tmp_yx1y[:, :, None, :], n_lu, axis=2)
            tmp_yx1y = np.sum(base["maps_YXvy"], axis=2)

    # Convert fraction of land to land area (m2)
    print('    Convert from "fraction of land" to "land area (m2)"')
    if do_harm:
        base_vegd_yx = np.sum(base["maps_YXv"][:, :, np.array(lu_names) != "BARREN"], axis=2)
        base_bare_yx = base["maps_YXv"][:, :, lu_names.index("BARREN")]
        if np.any(base_vegd_yx > 1 + 1e-15):
            raise RuntimeError("Floating-point error results in unacceptable base_vegd_YX!")
        base_bare_yx[base_vegd_yx > 1] = 0
        base["maps_YXv"][:, :, lu_names.index("BARREN")] = base_bare_yx
        base["maps_YXv"] = base["maps_YXv"] * land_area_yxv
    else:
        n_years_baseline = len(base["yearList"])
        base["maps_YXvy"] = base["maps_YXvy"] * np.repeat(
            land_area_yxv[:, :, :, None], n_years_baseline, axis=3
        )

    # Import base-year crop fractions
    print("    Import base_year crop fractions")
    if combine_crops:
        lpjg_crops = ["AllCrop"]
        n_crops_lpjg = len(lpjg_crops)
        base["varNames"][base["varNames"].index("CROPLAND")] = "AllCrop"
    else:
        print("      Read base_cropf")
        base_cropf = read_table_then2map(
            cfg.remap_cropf_file,
            verboseIfNoMat=False,
            force_mat_nosave=True,
            force_mat_save=False,
        )
        if "yearList" in base_cropf:
            print("      Get just base year")
            if do_harm:
                tmp = base_cropf["maps_YXvy"][:, :, :, base_cropf["yearList"] == cfg.base_year]
                base_cropf = {k: v for k, v in base_cropf.items() if k not in {"maps_YXvy", "yearList"}}
                base_cropf["maps_YXv"] = tmp.squeeze(axis=3)
                base_cropf["maps_YXv"][np.repeat(mask_yx[:, :, None], len(base_cropf["varNames"]), axis=2)] = np.nan
            else:
                n_years_baseline = len(base["yearList"])
                base_cropf["maps_YXvy"][
                    np.repeat(mask_yx[:, :, None, None], len(base_cropf["varNames"]), axis=2)
                ] = np.nan

        # Combine CC3G and CC4G into ExtraCrop
        if "CC3G" in base_cropf["varNames"] and "CC4G" in base_cropf["varNames"]:
            print("      Combine CC3G and CC4G into ExtraCrop")
            i_cc3g = base_cropf["varNames"].index("CC3G")
            i_cc4g = base_cropf["varNames"].index("CC4G")
            if "maps_YXv" in base_cropf:
                base_cropf["maps_YXv"][:, :, i_cc3g] = (
                    base_cropf["maps_YXv"][:, :, i_cc3g] + base_cropf["maps_YXv"][:, :, i_cc4g]
                )
                base_cropf["maps_YXv"] = np.delete(base_cropf["maps_YXv"], i_cc4g, axis=2)
            if "maps_YXvy" in base_cropf:
                base_cropf["maps_YXvy"][:, :, i_cc3g, :] = (
                    base_cropf["maps_YXvy"][:, :, i_cc3g, :]
                    + base_cropf["maps_YXvy"][:, :, i_cc4g, :]
                )
                base_cropf["maps_YXvy"] = np.delete(base_cropf["maps_YXvy"], i_cc4g, axis=2)
            base_cropf["varNames"].pop(i_cc4g)
            base_cropf["varNames"][i_cc3g] = "ExtraCrop"

        # Avoid sum(base_cropf) > 1
        if do_harm:
            print("      Avoid, to high precision, sum(base_cropf) > 1")
            tmp_yx = np.sum(base_cropf["maps_YXv"], axis=2)
            j = 0
            while np.any(tmp_yx - 1 > 3 * eps):
                j += 1
                if j > 50:
                    raise RuntimeError(
                        'Possible infinite loop in "Avoid, to high precision, sum(base_cropf) > 1"'
                    )
                base_cropf["maps_YXv"] = base_cropf["maps_YXv"] / np.repeat(
                    tmp_yx[:, :, None], base_cropf["maps_YXv"].shape[2], axis=2
                )
                tmp_yx = np.sum(base_cropf["maps_YXv"], axis=2)

        # Get previous crop types (remove irrigated duplicates)
        print("      Get previous crop types")
        tmp = list(base_cropf["varNames"])
        for name in list(base_cropf["varNames"]):
            if name.endswith("i") and name[:-1] in tmp:
                tmp.remove(name)
        lpjg_crops = tmp
        n_crops_lpjg = len(lpjg_crops)

        # Align yearList if not do_harm
        if not do_harm:
            n_years_baseline = len(base["yearList"])
            if "yearList" not in base_cropf:
                base_cropf["yearList"] = base["yearList"]
                base_cropf["maps_YXvy"] = np.repeat(
                    base_cropf["maps_YXv"][:, :, :, None], n_years_baseline, axis=3
                )
                base_cropf.pop("maps_YXv")
            else:
                year_is = [list(base_cropf["yearList"]).index(y) for y in base["yearList"]]
                base_cropf["maps_YXvy"] = base_cropf["maps_YXvy"][:, :, :, year_is]
                base_cropf["yearList"] = base["yearList"]

        # Irrigation intensity per crop
        print('      Get "irrigation intensity" as fraction of thisCrop that is irrigated')
        base_irrig: Dict[str, Any] = {"varNames": lpjg_crops}
        if do_harm:
            base_irrig["maps_YXv"] = np.full((*land_area_yx.shape, n_crops_lpjg), np.nan)
            for c, crop in enumerate(lpjg_crops):
                if crop == "ExtraCrop":
                    base_irrig["maps_YXv"][:, :, c] = 0
                    continue
                crop_i = crop + "i"
                crop_r = base_cropf["maps_YXv"][:, :, base_cropf["varNames"].index(crop)]
                crop_i_yx = base_cropf["maps_YXv"][:, :, base_cropf["varNames"].index(crop_i)]
                crop_yx = crop_r + crop_i_yx
                tmp_yx = crop_i_yx / crop_yx
                tmp_yx[crop_yx == 0] = np.nan
                base_irrig["maps_YXv"][:, :, c] = tmp_yx
        else:
            base_irrig["maps_YXvy"] = np.full(
                (*land_area_yx.shape, n_crops_lpjg, len(base["yearList"])), np.nan
            )
            for c, crop in enumerate(lpjg_crops):
                if crop == "ExtraCrop":
                    base_irrig["maps_YXvy"][:, :, c, :] = 0
                    continue
                crop_i = crop + "i"
                crop_r = base_cropf["maps_YXvy"][:, :, base_cropf["varNames"].index(crop), :]
                crop_i_yx = base_cropf["maps_YXvy"][:, :, base_cropf["varNames"].index(crop_i), :]
                crop_yx = crop_r + crop_i_yx
                tmp_yx = crop_i_yx / crop_yx
                tmp_yx[crop_yx == 0] = np.nan
                base_irrig["maps_YXvy"][:, :, c, :] = tmp_yx

        # Combine irrigated and rainfed
        print("      Combine irrigated and rainfed")
        base_cropf = _combine_irrig(base_cropf, lpjg_crops)

        # Read baseline fertilization
        print("      Read baseline fertilization")
        base_nfert = read_table_then2map(
            cfg.remap_nfert_file,
            verboseIfNoMat=False,
            force_mat_nosave=True,
            force_mat_save=False,
        )
        if do_harm and "yearList" in base_nfert:
            tmp = base_nfert["maps_YXvy"][:, :, :, base_nfert["yearList"] == cfg.base_year]
            base_nfert = {k: v for k, v in base_nfert.items() if k not in {"maps_YXvy", "yearList"}}
            base_nfert["maps_YXv"] = tmp.squeeze(axis=3)
        elif not do_harm:
            n_years_baseline = len(base["yearList"])
            if "yearList" not in base_nfert:
                base_nfert["yearList"] = base["yearList"]
                base_nfert["maps_YXvy"] = np.repeat(
                    base_nfert["maps_YXv"][:, :, :, None], n_years_baseline, axis=3
                )
                base_nfert.pop("maps_YXv")
            else:
                year_is = [list(base_nfert["yearList"]).index(y) for y in base["yearList"]]
                base_nfert["maps_YXvy"] = base_nfert["maps_YXvy"][:, :, :, year_is]
                base_nfert["yearList"] = base["yearList"]

        # Add ExtraCrop if needed
        if "ExtraCrop" in lpjg_crops and "ExtraCrop" not in base_nfert["varNames"]:
            print("      Add ExtraCrop")
            base_nfert["varNames"].append("ExtraCrop")
            if do_harm:
                base_nfert["maps_YXv"] = np.concatenate(
                    [base_nfert["maps_YXv"], np.zeros((*land_area_yx.shape, 1))], axis=2
                )
            else:
                base_nfert["maps_YXvy"] = np.concatenate(
                    [
                        base_nfert["maps_YXvy"],
                        np.zeros((*land_area_yx.shape, 1, len(base["yearList"]))),
                    ],
                    axis=2,
                )

        # Reorder nfert to LPJ crops
        idx = [base_nfert["varNames"].index(c) for c in lpjg_crops]
        if do_harm:
            base_nfert["maps_YXv"] = base_nfert["maps_YXv"][:, :, idx]
        else:
            base_nfert["maps_YXvy"] = base_nfert["maps_YXvy"][:, :, idx, :]
        base_nfert["varNames"] = list(lpjg_crops)

        # Convert base to crop-level instead of CROPLAND
        base_orig = base
        base_cropa = base_cropf
        if do_harm:
            print("      Convert base to have individual crops instead of CROPLAND")
            i_crop = base_orig["varNames"].index("CROPLAND")
            base_cropa["maps_YXv"] = base_cropf["maps_YXv"] * base_orig["maps_YXv"][:, :, i_crop][:, :, None]
            base = {
                "varNames": base_cropa["varNames"] + [
                    v for v in base_orig["varNames"] if v != "CROPLAND"
                ],
                "maps_YXv": np.concatenate(
                    [base_cropa["maps_YXv"], base_orig["maps_YXv"][:, :, i_crop + 1 :]], axis=2
                ),
            }
        else:
            i_crop = base_orig["varNames"].index("CROPLAND")
            base_cropa["maps_YXvy"] = base_cropf["maps_YXvy"] * base_orig["maps_YXvy"][:, :, i_crop, :][
                :, :, None, :
            ]
            base = {
                "varNames": base_cropa["varNames"] + [
                    v for v in base_orig["varNames"] if v != "CROPLAND"
                ],
                "maps_YXvy": np.concatenate(
                    [base_cropa["maps_YXvy"], base_orig["maps_YXvy"][:, :, i_crop + 1 :, :]], axis=2
                ),
            }

    # Info
    (
        lu_names,
        max_len_lu,
        n_lu,
        not_bare,
        is_agri,
        lu_names_agri,
        is_agri_is_past,
        is_crop,
    ) = _get_varname_info(base["varNames"])

    # Convert NaNs to zeros for addition
    print("    Convert NaNs to zeros for addition")
    if do_harm:
        base["maps_YXv"][np.isnan(base["maps_YXv"])] = 0
        if not combine_crops:
            base_nfert["maps_YXv"][base["maps_YXv"][:, :, is_crop] == 0] = 0
            base_irrig["maps_YXv"][base["maps_YXv"][:, :, is_crop] == 0] = 0
    else:
        base["maps_YXvy"][np.isnan(base["maps_YXvy"])] = 0
        if not combine_crops:
            base_nfert["maps_YXvy"][base["maps_YXvy"][:, :, is_crop, :] == 0] = 0
            base_irrig["maps_YXvy"][base["maps_YXvy"][:, :, is_crop, :] == 0] = 0

    # Baseline veg/bare fractions
    base_vegd_yx, base_bare_yx, base_vegd_frac_yx, base_bare_frac_yx = _get_bl_frac_vegd_bare(
        base, land_area_yx, not_bare, lu_names, do_harm
    )

    # Aggregate to 2-degree
    print("    Aggregate from 0.5-degree to 2-degree")
    base_2deg = {"varNames": base["varNames"]}
    base_2deg_nfert = {"varNames": lpjg_crops} if not combine_crops else None
    base_2deg_irrig = {"varNames": lpjg_crops} if not combine_crops else None
    if do_harm:
        base_2deg["maps_YXv"] = aggregate(base["maps_YXv"], 0.5, 2)
        if not combine_crops:
            base_2deg_nfert["maps_YXv"] = _aggregate_mgmt(
                base_nfert["maps_YXv"], base["maps_YXv"][:, :, is_crop], 0.5, 2
            )
            base_2deg_irrig["maps_YXv"] = _aggregate_mgmt(
                base_irrig["maps_YXv"], base["maps_YXv"][:, :, is_crop], 0.5, 2
            )
    else:
        base_2deg["maps_YXvy"] = aggregate(base["maps_YXvy"], 0.5, 2)
        if not combine_crops:
            base_2deg_nfert["maps_YXvy"] = _aggregate_mgmt(
                base_nfert["maps_YXvy"], base["maps_YXvy"][:, :, is_crop, :], 0.5, 2
            )
            base_2deg_irrig["maps_YXvy"] = _aggregate_mgmt(
                base_irrig["maps_YXvy"], base["maps_YXvy"][:, :, is_crop, :], 0.5, 2
            )
    land_area_2deg_yx = aggregate(land_area_yx, 0.5, 2)

    # Validate mgmt inputs
    print("    Do not allow invalid management inputs")
    if not combine_crops:
        if do_harm:
            if np.any(base_2deg_nfert["maps_YXv"] < 0):
                raise RuntimeError("Negative value(s) in base_2deg_nfert.maps_YXv!")
            if np.any(base_2deg_irrig["maps_YXv"] < 0):
                raise RuntimeError("Negative value(s) in base_2deg_irrig.maps_YXv!")
            if np.any(base_2deg_irrig["maps_YXv"] > 1 + 1e-6):
                raise RuntimeError("Value(s) >1 in base_2deg_irrig.maps_YXv!")
        else:
            if np.any(base_2deg_nfert["maps_YXvy"] < 0):
                raise RuntimeError("Negative value(s) in base_2deg_nfert.maps_YXvy!")
            if np.any(base_2deg_irrig["maps_YXvy"] < 0):
                raise RuntimeError("Negative value(s) in base_2deg_irrig.maps_YXvy!")
            if np.any(base_2deg_irrig["maps_YXvy"] > 1 + 1e-6):
                raise RuntimeError("Value(s) >1 in base_2deg_irrig.maps_YXvy!")

    # Interpolate management inputs (doHarm only)
    print("    Interpolate management inputs")
    if not combine_crops and do_harm:
        base_2deg_nfert["maps_YXv"] = _interpolate_mgmt(
            base_2deg_nfert["maps_YXv"],
            base["maps_YXv"][:, :, is_crop],
            land_area_2deg_yx,
            lpjg_crops,
            cfg.inpaint_method,
        )
        base_2deg_irrig["maps_YXv"] = _interpolate_mgmt(
            base_2deg_irrig["maps_YXv"],
            base["maps_YXv"][:, :, is_crop],
            land_area_2deg_yx,
            lpjg_crops,
            cfg.inpaint_method,
        )

    # Baseline LU 2-deg veg/bare fractions
    print("    Get baseline LU 2-deg fraction that is vegetated, barren")
    if do_harm:
        base_2deg_vegd_yx = np.sum(base_2deg["maps_YXv"][:, :, not_bare], axis=2)
        base_2deg_bare_yx = base_2deg["maps_YXv"][:, :, lu_names.index("BARREN")]
    else:
        base_2deg_vegd_yx = np.sum(base_2deg["maps_YXvy"][:, :, not_bare, 0], axis=2)
        base_2deg_bare_yx = base_2deg["maps_YXvy"][:, :, lu_names.index("BARREN"), 0]
    base_2deg_vegd_frac_yx = base_2deg_vegd_yx / land_area_2deg_yx
    base_2deg_bare_frac_yx = base_2deg_bare_yx / land_area_2deg_yx

    # Other info
    print("    Get other info")
    base_nfert_tot = None
    base_nfert_tot_2deg = None
    base_irrig_tot = None
    base_irrig_tot_2deg = None
    if not combine_crops:
        if do_harm:
            base_nfert_tot_2deg = {"maps_YXv": base_2deg_nfert["maps_YXv"] * base_2deg["maps_YXv"][:, :, is_crop]}
            base_irrig_tot_2deg = {"maps_YXv": base_2deg_irrig["maps_YXv"] * base_2deg["maps_YXv"][:, :, is_crop]}
        else:
            base_nfert_tot = {"maps_YXvy": base_nfert["maps_YXvy"] * base["maps_YXvy"][:, :, is_crop, :]}
            base_nfert_tot_2deg = {
                "maps_YXvy": base_2deg_nfert["maps_YXvy"] * base_2deg["maps_YXvy"][:, :, is_crop, :]
            }
            base_irrig_tot = {"maps_YXvy": base_irrig["maps_YXvy"] * base["maps_YXvy"][:, :, is_crop, :]}
            base_irrig_tot_2deg = {
                "maps_YXvy": base_2deg_irrig["maps_YXvy"] * base_2deg["maps_YXvy"][:, :, is_crop, :]
            }

    # lats/lons and list2map
    list2map_2deg = np.flatnonzero(land_area_2deg_yx.ravel(order="F") > 0) + 1
    lons_2deg = lons_map_2deg.ravel(order="F")[list2map_2deg - 1]
    lats_2deg = lats_map_2deg.ravel(order="F")[list2map_2deg - 1]
    list2map = np.flatnonzero(land_area_yx.ravel(order="F") > 0) + 1
    lons = lons_map.ravel(order="F")[list2map - 1]
    lats = lats_map.ravel(order="F")[list2map - 1]

    # PLUM crop types (from LandUse.txt header)
    file_in = os.path.join(dir1_base, "LandUse.txt")
    T = read_table(file_in, verboseIfNoMat=False, dont_save_MAT=True, do_save_MAT=False)
    plum_crops = [
        c
        for c in T.columns
        if ("_A" in c)
        and ("ruminants" not in c)
        and ("monogastrics" not in c)
        and ("pasture" not in c)
    ]
    plum_crops = [c.replace("_A", "") for c in plum_crops]

    # PLUM to LPJG mapping
    if combine_crops:
        plum_to_lpjg = []
    else:
        plum_to_lpjg = [None] * len(lpjg_crops)
        _set_map(plum_to_lpjg, lpjg_crops, "CerealsC3", "wheat")
        _set_map(plum_to_lpjg, lpjg_crops, "CerealsC4", "maize")
        _set_map(plum_to_lpjg, lpjg_crops, "Rice", "rice")
        if "Oilcrops" in lpjg_crops:
            _set_map(plum_to_lpjg, lpjg_crops, "Oilcrops", "oilcrops")
        else:
            has_split = "OilNfix" in lpjg_crops and "OilOther" in lpjg_crops
            if not has_split:
                raise RuntimeError("Oilcrops not in LPJcrops and split OilNfix/OilOther not found.")
            _set_map(plum_to_lpjg, lpjg_crops, "OilNfix", "oilcropsNFix")
            _set_map(plum_to_lpjg, lpjg_crops, "OilOther", "oilcropsOther")
        _set_map(plum_to_lpjg, lpjg_crops, "Pulses", "pulses")
        _set_map(plum_to_lpjg, lpjg_crops, "StarchyRoots", "starchyRoots")
        if "FruitVeg" in lpjg_crops:
            _set_map(plum_to_lpjg, lpjg_crops, "FruitVeg", "fruitveg")
        if "FruitAndVeg" in lpjg_crops:
            _set_map(plum_to_lpjg, lpjg_crops, "FruitAndVeg", "fruitveg")
        if "Sugar" in lpjg_crops:
            _set_map(plum_to_lpjg, lpjg_crops, "Sugar", "sugar")
        _set_map(plum_to_lpjg, lpjg_crops, "Miscanthus", "energycrops")
        _set_map(plum_to_lpjg, lpjg_crops, "ExtraCrop", "setaside")
        _check_plum_map(plum_to_lpjg, plum_crops, cfg.fruitveg_sugar_2oil)

    # No-unveg check
    if not cfg.allow_unveg:
        _check_no_unveg(base, land_area_yx, "Half-deg")
        _check_no_unveg(base_2deg, land_area_2deg_yx, "2-deg")

    print("Done importing reference data.")

    return {
        "lons_map": lons_map,
        "lats_map": lats_map,
        "lons_map_2deg": lons_map_2deg,
        "lats_map_2deg": lats_map_2deg,
        "landArea_YX": land_area_yx,
        "landArea_2deg_YX": land_area_2deg_yx,
        "gcelArea_YX": gcel_area_yx,
        "resFrac_terr_YX": res_frac_terr_yx,
        "resFrac_prot_YX": res_frac_prot_yx,
        "base": base,
        "base_2deg": base_2deg,
        "base_nfert": base_nfert if not combine_crops else None,
        "base_irrig": base_irrig if not combine_crops else None,
        "base_2deg_nfert": base_2deg_nfert if not combine_crops else None,
        "base_2deg_irrig": base_2deg_irrig if not combine_crops else None,
        "base_nfert_tot": base_nfert_tot,
        "base_nfert_tot_2deg": base_nfert_tot_2deg,
        "base_irrig_tot": base_irrig_tot,
        "base_irrig_tot_2deg": base_irrig_tot_2deg,
        "base_vegd_YX": base_vegd_yx,
        "base_bare_YX": base_bare_yx,
        "base_vegdFrac_YX": base_vegd_frac_yx,
        "base_bareFrac_YX": base_bare_frac_yx,
        "base_2deg_vegd_YX": base_2deg_vegd_yx,
        "base_2deg_bare_YX": base_2deg_bare_yx,
        "base_2deg_vegdFrac_YX": base_2deg_vegd_frac_yx,
        "base_2deg_bareFrac_YX": base_2deg_bare_frac_yx,
        "lu_names": lu_names,
        "lpjg_crops": lpjg_crops if not combine_crops else ["AllCrop"],
        "plum_crops": plum_crops,
        "plum_to_lpjg": plum_to_lpjg,
        "list2map": list2map,
        "list2map_2deg": list2map_2deg,
        "lons": lons,
        "lats": lats,
        "lons_2deg": lons_2deg,
        "lats_2deg": lats_2deg,
        "mask_YX": mask_yx,
        "notBare": not_bare,
        "isAgri": is_agri,
        "isCrop": is_crop,
        "maxLength_LUnames": max_len_lu,
    }


def _get_varname_info(
    var_names: Sequence[str],
) -> Tuple[Sequence[str], int, int, np.ndarray, np.ndarray, Sequence[str], np.ndarray, np.ndarray]:
    lu_names = list(var_names)
    max_len_lu = max(len(v) for v in lu_names)
    n_lu = len(lu_names)
    not_bare = np.array([v != "BARREN" for v in lu_names])
    is_agri = np.array([v not in {"NATURAL", "BARREN"} for v in lu_names])
    lu_names_agri = [v for v, is_a in zip(lu_names, is_agri) if is_a]
    is_agri_is_past = np.array([v == "PASTURE" for v in lu_names_agri])
    is_crop = np.array([v not in {"NATURAL", "PASTURE", "BARREN"} for v in lu_names])
    return lu_names, max_len_lu, n_lu, not_bare, is_agri, lu_names_agri, is_agri_is_past, is_crop


def _get_bl_frac_vegd_bare(base, land_area_yx, not_bare, lu_names, do_harm):
    print("    Get baseline fraction that is vegetated, barren")
    if do_harm:
        base_vegd_yx = np.sum(base["maps_YXv"][:, :, not_bare], axis=2)
        base_bare_yx = base["maps_YXv"][:, :, lu_names.index("BARREN")]
    else:
        base_vegd_yx = np.sum(base["maps_YXvy"][:, :, not_bare, 0], axis=2)
        base_bare_yx = base["maps_YXvy"][:, :, lu_names.index("BARREN"), 0]
    base_vegd_frac_yx = base_vegd_yx / land_area_yx
    base_bare_frac_yx = base_bare_yx / land_area_yx
    return base_vegd_yx, base_bare_yx, base_vegd_frac_yx, base_bare_frac_yx


def _combine_irrig(cropf_in, lpjg_crops):
    n_crops_lpjg = len(lpjg_crops)
    any_irrigated = cropf_in["varNames"] != lpjg_crops
    if not any_irrigated:
        return cropf_in
    cropf_out = {"varNames": list(lpjg_crops)}
    if "maps_YXv" in cropf_in:
        cropf_out["maps_YXv"] = np.full(
            (cropf_in["maps_YXv"].shape[0], cropf_in["maps_YXv"].shape[1], n_crops_lpjg), np.nan
        )
        for c, crop in enumerate(lpjg_crops):
            mask = [crop in v for v in cropf_in["varNames"]]
            cropf_out["maps_YXv"][:, :, c] = np.sum(cropf_in["maps_YXv"][:, :, mask], axis=2)
    else:
        cropf_out["yearList"] = cropf_in["yearList"]
        cropf_out["maps_YXvy"] = np.full(
            (
                cropf_in["maps_YXvy"].shape[0],
                cropf_in["maps_YXvy"].shape[1],
                n_crops_lpjg,
                cropf_in["maps_YXvy"].shape[3],
            ),
            np.nan,
        )
        for c, crop in enumerate(lpjg_crops):
            mask = [crop in v for v in cropf_in["varNames"]]
            cropf_out["maps_YXvy"][:, :, c, :] = np.sum(cropf_in["maps_YXvy"][:, :, mask, :], axis=2)
    return cropf_out


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


def _check_no_unveg(base, land_area_yx, label):
    var_names = list(base.get("varNames", []))
    not_bare = [i for i, v in enumerate(var_names) if v != "BARREN"]
    if "maps_YXv" in base:
        bad = np.sum(base["maps_YXv"][:, :, not_bare], axis=2) == 0
    else:
        bad = np.any(
            np.sum(base["maps_YXvy"][:, :, not_bare, :], axis=2) == 0, axis=2
        )
    n_bad = int(np.sum(bad & (land_area_yx > 0)))
    if n_bad > 0:
        raise RuntimeError(
            f"{label} baseline has {n_bad} non-vegetated gridcells where landArea>0!"
        )


def _set_map(plum_to_lpjg, lpjg_crops, lpjg_name, plum_name):
    if lpjg_name in lpjg_crops:
        plum_to_lpjg[lpjg_crops.index(lpjg_name)] = plum_name


def _check_plum_map(plum_to_lpjg, plum_crops, fruitveg_sugar_2oil):
    mapped = set(v for v in plum_to_lpjg if v)
    missing = [c for c in plum_crops if c not in mapped]
    if missing:
        raise RuntimeError(f"Unmapped PLUM crops: {missing} (fruitveg_sugar_2oil={fruitveg_sugar_2oil})")
