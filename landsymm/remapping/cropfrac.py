"""Crop fraction processing (stub).

Pseudocode:
1) Read MIRCA harvested area tables (RF/IR).
2) Build combined crop list + mappings from get_remapv2_keys(this_ver).
3) Aggregate MIRCA crops into output crop list (RF/IR).
4) Apply force_all_rainfed if enabled.
5) Compute ExtraCrop from ignored crops + setAside.
6) Normalize crop fractions to sum to 1.
"""

# TODO:
# - implement MIRCA read and mapping parity
# - preserve crop ordering and RF/IR handling
from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

from landsymm.common.lpjg_io import read_table_then2map
from landsymm.common.mapping_tools import yxz_to_xz


def get_remapv2_keys(this_ver: str):
    """Return crop mapping for this_ver (MATLAB get_remapv2_keys.m)."""
    out_names = [
        "20180105b",
        "20180206",
        "20180210",
        "20180212",
        "20180214",
        "20190216",
        "20180301ani",
        "WithFruitVegSugar_a",
        "WithFruitVegSugar_b",
        "WithFruitVegSugar_b_2oil",
        "WithFruitVeg_sepSugar",
        "WithFruitVeg_sepSugar_sepOil",
        "WithFruitVeg_sepOil_combSugar",
        "WithFruitVeg_sepSugar_sepOilInclPalm",
        "WithFruitVeg_sepSugar_sepOil_sepC3",
        "ggcmi5",
        "ggcmi5_preBNF",
        "jianyong01",
        "jianyong01b",
        "ani01",
        "20200928",
    ]
    if this_ver not in out_names:
        raise RuntimeError(f"this_ver {this_ver} not recognized!")

    out_keys: dict[str, list[list[str]]] = {}
    out_list: dict[str, list[str]] = {}
    out_ignores: dict[str, list[str]] = {}

    # 20180105b
    this = "20180105b"
    crops = ["TeWW", "TeSW", "TeCo", "TrRi"]
    out_keys[this] = [
        ["Wheat", "Barley"],
        ["Pulses", "Potatoes", "Sugarbeet", "Cassava", "Sunflower", "Soybeans", "GroundnutsPeanuts", "RapeseedCanola"],
        ["Maize", "Millet", "Sorghum"],
        ["Rice"],
    ]
    out_list[this] = crops
    out_ignores[this] = [
        "Sugarcane",
        "Oilpalm",
        "Citrus",
        "Datepalm",
        "GrapesVine",
        "Cotton",
        "Cocoa",
        "Coffee",
        "OtherAnnuals",
        "OtherPerennials",
        "FodderGrasses",
        "Rye",
    ]

    # 20180206
    this = "20180206"
    crops = ["Wheat", "Maize", "Sorghum", "Rice", "Soybeans", "Pulses"]
    out_keys[this] = [["Wheat"], ["Maize"], ["Sorghum"], ["Rice"], ["Soybeans"], ["Pulses"]]
    out_list[this] = crops
    out_ignores[this] = [
        "Sugarcane",
        "Oilpalm",
        "Citrus",
        "Datepalm",
        "GrapesVine",
        "Cotton",
        "Cocoa",
        "Coffee",
        "OtherAnnuals",
        "OtherPerennials",
        "FodderGrasses",
        "Rye",
        "Barley",
        "Millet",
        "Sunflower",
        "GroundnutsPeanuts",
        "RapeseedCanola",
        "Potatoes",
        "Sugarbeet",
        "Cassava",
    ]

    # 20180210
    this = "20180210"
    crops = ["CerealsC3", "CerealsC4", "Rice", "Oilcrops", "Pulses", "StarchyRoots"]
    out_keys[this] = [
        ["Wheat", "Barley"],
        ["Maize", "Millet", "Sorghum"],
        ["Rice"],
        ["Sunflower", "Soybeans", "GroundnutsPeanuts", "RapeseedCanola"],
        ["Pulses"],
        ["Potatoes", "Sugarbeet", "Cassava"],
    ]
    out_list[this] = crops
    out_ignores[this] = [
        "Sugarcane",
        "Oilpalm",
        "Citrus",
        "Datepalm",
        "GrapesVine",
        "Cotton",
        "Cocoa",
        "Coffee",
        "OtherAnnuals",
        "OtherPerennials",
        "FodderGrasses",
        "Rye",
    ]

    # 20180212
    this = "20180212"
    crops = ["Wheat", "Maize", "Sorghum", "Rice"]
    out_keys[this] = [["Wheat"], ["Maize"], ["Sorghum"], ["Rice"]]
    out_list[this] = crops
    out_ignores[this] = [
        "Sugarcane",
        "Oilpalm",
        "Citrus",
        "Datepalm",
        "GrapesVine",
        "Cotton",
        "Cocoa",
        "Coffee",
        "OtherAnnuals",
        "OtherPerennials",
        "FodderGrasses",
        "Rye",
        "Barley",
        "Millet",
        "Sunflower",
        "GroundnutsPeanuts",
        "RapeseedCanola",
        "Potatoes",
        "Sugarbeet",
        "Cassava",
        "Soybeans",
        "Pulses",
    ]

    # 20180214
    this = "20180214"
    crops = ["CerealsC3", "CerealsC4", "Rice", "Oilcrops", "Pulses", "StarchyRoots"]
    out_keys[this] = [
        ["Wheat", "Barley", "Rye"],
        ["Maize", "Millet", "Sorghum"],
        ["Rice"],
        ["Sunflower", "Soybeans", "GroundnutsPeanuts", "RapeseedCanola", "Oilpalm"],
        ["Pulses"],
        ["Potatoes", "Sugarbeet", "Cassava"],
    ]
    out_list[this] = crops
    out_ignores[this] = [
        "Sugarcane",
        "Citrus",
        "Datepalm",
        "GrapesVine",
        "Cotton",
        "Cocoa",
        "Coffee",
        "OtherAnnuals",
        "OtherPerennials",
        "FodderGrasses",
    ]

    # 20190216
    this = "20190216"
    crops = ["CerealsC3", "CerealsC4", "Rice", "Oilcrops", "Pulses", "StarchyRoots"]
    out_keys[this] = [
        ["Wheat", "Barley", "Rye"],
        ["Maize", "Millet", "Sorghum"],
        ["Rice"],
        ["Sunflower", "Soybeans", "GroundnutsPeanuts", "RapeseedCanola", "Oilpalm"],
        ["Pulses"],
        ["Potatoes", "Cassava"],
    ]
    out_list[this] = crops
    out_ignores[this] = [
        "Sugarcane",
        "Citrus",
        "Datepalm",
        "GrapesVine",
        "Cotton",
        "Cocoa",
        "Coffee",
        "OtherAnnuals",
        "OtherPerennials",
        "FodderGrasses",
        "Sugarbeet",
    ]

    # 20180301ani
    this = "20180301ani"
    crops = ["TeWW", "TeSW", "TeCo", "TrRi"]
    out_keys[this] = [
        ["Wheat", "Barley", "Rye"],
        ["Sunflower", "Soybeans", "GroundnutsPeanuts", "RapeseedCanola", "Oilpalm", "Pulses", "Potatoes", "Sugarbeet", "Cassava"],
        ["Maize", "Millet", "Sorghum"],
        ["Rice"],
    ]
    out_list[this] = crops
    out_ignores[this] = [
        "Sugarcane",
        "Citrus",
        "Datepalm",
        "GrapesVine",
        "Cotton",
        "Cocoa",
        "Coffee",
        "OtherAnnuals",
        "OtherPerennials",
        "FodderGrasses",
    ]

    # WithFruitVegSugar_a
    this = "WithFruitVegSugar_a"
    crops = ["CerealsC3", "CerealsC4", "Rice", "Oilcrops", "Pulses", "StarchyRoots", "DateCitGrape", "Sugar"]
    out_keys[this] = [
        ["Wheat", "Barley", "Rye"],
        ["Maize", "Millet", "Sorghum"],
        ["Rice"],
        ["Sunflower", "Soybeans", "GroundnutsPeanuts", "RapeseedCanola", "Oilpalm"],
        ["Pulses"],
        ["Potatoes", "Cassava"],
        ["Datepalm", "Citrus", "GrapesVine"],
        ["Sugarbeet", "Sugarcane"],
    ]
    out_list[this] = crops
    out_ignores[this] = ["Cotton", "Cocoa", "Coffee", "FodderGrasses", "OtherAnnuals", "OtherPerennials"]

    # WithFruitVegSugar_b
    this = "WithFruitVegSugar_b"
    crops = ["CerealsC3", "CerealsC4", "Rice", "Oilcrops", "Pulses", "StarchyRoots", "FruitAndVeg", "Sugar"]
    out_keys[this] = [
        ["Wheat", "Barley", "Rye"],
        ["Maize", "Millet", "Sorghum"],
        ["Rice"],
        ["Sunflower", "Soybeans", "GroundnutsPeanuts", "RapeseedCanola", "Oilpalm"],
        ["Pulses"],
        ["Potatoes", "Cassava"],
        ["Datepalm", "Citrus", "GrapesVine", "OtherAnnuals", "OtherPerennials"],
        ["Sugarbeet", "Sugarcane"],
    ]
    out_list[this] = crops
    out_ignores[this] = ["Cotton", "Cocoa", "Coffee", "FodderGrasses"]

    # WithFruitVegSugar_b_2oil
    this = "WithFruitVegSugar_b_2oil"
    crops = ["CerealsC3", "CerealsC4", "Rice", "Oilcrops", "Pulses", "StarchyRoots"]
    out_keys[this] = [
        ["Wheat", "Barley", "Rye"],
        ["Maize", "Millet", "Sorghum"],
        ["Rice"],
        ["Sunflower", "Soybeans", "GroundnutsPeanuts", "RapeseedCanola", "Oilpalm", "Datepalm", "Citrus", "GrapesVine", "OtherAnnuals", "OtherPerennials", "Sugarbeet", "Sugarcane"],
        ["Pulses"],
        ["Potatoes", "Cassava"],
    ]
    out_list[this] = crops
    out_ignores[this] = ["Cotton", "Cocoa", "Coffee", "FodderGrasses"]

    # WithFruitVeg_sepSugar
    this = "WithFruitVeg_sepSugar"
    crops = ["CerealsC3", "CerealsC4", "Rice", "Oilcrops", "Pulses", "StarchyRoots", "FruitAndVeg", "Sugarbeet", "Sugarcane"]
    out_keys[this] = [
        ["Wheat", "Barley", "Rye"],
        ["Maize", "Millet", "Sorghum"],
        ["Rice"],
        ["Sunflower", "Soybeans", "GroundnutsPeanuts", "RapeseedCanola", "Oilpalm"],
        ["Pulses"],
        ["Potatoes", "Cassava"],
        ["Datepalm", "Citrus", "GrapesVine", "OtherAnnuals", "OtherPerennials"],
        ["Sugarbeet"],
        ["Sugarcane"],
    ]
    out_list[this] = crops
    out_ignores[this] = ["Cotton", "Cocoa", "Coffee", "FodderGrasses"]

    # WithFruitVeg_sepSugar_sepOil
    this = "WithFruitVeg_sepSugar_sepOil"
    crops = ["CerealsC3", "CerealsC4", "Rice", "OilNfix", "OilOther", "Pulses", "StarchyRoots", "FruitAndVeg", "Sugarbeet", "Sugarcane"]
    out_keys[this] = [
        ["Wheat", "Barley", "Rye"],
        ["Maize", "Millet", "Sorghum"],
        ["Rice"],
        ["Soybeans", "GroundnutsPeanuts"],
        ["Sunflower", "RapeseedCanola", "Oilpalm"],
        ["Pulses"],
        ["Potatoes", "Cassava"],
        ["Datepalm", "Citrus", "GrapesVine", "OtherAnnuals", "OtherPerennials"],
        ["Sugarbeet"],
        ["Sugarcane"],
    ]
    out_list[this] = crops
    out_ignores[this] = ["Cotton", "Cocoa", "Coffee", "FodderGrasses"]

    # WithFruitVeg_sepOil_combSugar
    this = "WithFruitVeg_sepOil_combSugar"
    crops = ["CerealsC3", "CerealsC4", "Rice", "OilNfix", "OilOther", "Pulses", "StarchyRoots", "FruitAndVeg", "Sugar"]
    out_keys[this] = [
        ["Wheat", "Barley", "Rye"],
        ["Maize", "Millet", "Sorghum"],
        ["Rice"],
        ["Soybeans", "GroundnutsPeanuts"],
        ["Sunflower", "RapeseedCanola", "Oilpalm"],
        ["Pulses"],
        ["Potatoes", "Cassava"],
        ["Datepalm", "Citrus", "GrapesVine", "OtherAnnuals", "OtherPerennials"],
        ["Sugarbeet", "Sugarcane"],
    ]
    out_list[this] = crops
    out_ignores[this] = ["Cotton", "Cocoa", "Coffee", "FodderGrasses"]

    # WithFruitVeg_sepSugar_sepOilInclPalm
    this = "WithFruitVeg_sepSugar_sepOilInclPalm"
    crops = ["CerealsC3", "CerealsC4", "Rice", "OilNfix", "OilOther", "OilPalm", "Pulses", "StarchyRoots", "FruitAndVeg", "Sugarbeet", "Sugarcane"]
    out_keys[this] = [
        ["Wheat", "Barley", "Rye"],
        ["Maize", "Millet", "Sorghum"],
        ["Rice"],
        ["Soybeans", "GroundnutsPeanuts"],
        ["Sunflower", "RapeseedCanola"],
        ["Oilpalm"],
        ["Pulses"],
        ["Potatoes", "Cassava"],
        ["Datepalm", "Citrus", "GrapesVine", "OtherAnnuals", "OtherPerennials"],
        ["Sugarbeet"],
        ["Sugarcane"],
    ]
    out_list[this] = crops
    out_ignores[this] = ["Cotton", "Cocoa", "Coffee", "FodderGrasses"]

    # WithFruitVeg_sepSugar_sepOil_sepC3
    this = "WithFruitVeg_sepSugar_sepOil_sepC3"
    crops = ["CerealsC3s", "CerealsC3w", "CerealsC4", "Rice", "OilNfix", "OilOther", "Pulses", "StarchyRoots", "FruitAndVeg", "Sugarbeet", "Sugarcane"]
    out_keys[this] = [
        ["Wheat", "Barley", "Rye"],
        ["Wheat", "Barley", "Rye"],
        ["Maize", "Millet", "Sorghum"],
        ["Rice"],
        ["Soybeans", "GroundnutsPeanuts"],
        ["Sunflower", "RapeseedCanola", "Oilpalm"],
        ["Pulses"],
        ["Potatoes", "Cassava"],
        ["Datepalm", "Citrus", "GrapesVine", "OtherAnnuals", "OtherPerennials"],
        ["Sugarbeet"],
        ["Sugarcane"],
    ]
    out_list[this] = crops
    out_ignores[this] = ["Cotton", "Cocoa", "Coffee", "FodderGrasses"]

    # ggcmi5
    this = "ggcmi5"
    crops = ["CerealsC3s", "CerealsC3w", "CerealsC4", "Rice", "Oilcrops"]
    out_keys[this] = [["Wheat"], ["Wheat"], ["Maize"], ["Rice"], ["Soybeans"]]
    out_list[this] = crops
    out_ignores[this] = [
        "Barley",
        "Rye",
        "Millet",
        "Sorghum",
        "GroundnutsPeanuts",
        "Sunflower",
        "RapeseedCanola",
        "Oilpalm",
        "Pulses",
        "Potatoes",
        "Cassava",
        "Datepalm",
        "Citrus",
        "GrapesVine",
        "OtherAnnuals",
        "OtherPerennials",
        "Sugarbeet",
        "Sugarcane",
        "Cotton",
        "Cocoa",
        "Coffee",
        "FodderGrasses",
    ]

    # ggcmi5_preBNF
    this = "ggcmi5_preBNF"
    crops = ["CerealsC3s", "CerealsC3w", "CerealsC4", "Rice"]
    out_keys[this] = [["Wheat"], ["Wheat"], ["Maize"], ["Rice"]]
    out_list[this] = crops
    out_ignores[this] = [
        "Soybeans",
        "Barley",
        "Rye",
        "Millet",
        "Sorghum",
        "GroundnutsPeanuts",
        "Sunflower",
        "RapeseedCanola",
        "Oilpalm",
        "Pulses",
        "Potatoes",
        "Cassava",
        "Datepalm",
        "Citrus",
        "GrapesVine",
        "OtherAnnuals",
        "OtherPerennials",
        "Sugarbeet",
        "Sugarcane",
        "Cotton",
        "Cocoa",
        "Coffee",
        "FodderGrasses",
    ]

    # jianyong01
    this = "jianyong01"
    crops = ["Wheat", "Maize", "Sorghum", "Rice", "Soybean", "FabaBean"]
    out_keys[this] = [["Wheat"], ["Maize"], ["Sorghum"], ["Rice"], ["Soybeans"], ["Pulses"]]
    out_list[this] = crops
    out_ignores[this] = [
        "Sugarcane",
        "Oilpalm",
        "Citrus",
        "Datepalm",
        "GrapesVine",
        "Cotton",
        "Cocoa",
        "Coffee",
        "OtherAnnuals",
        "OtherPerennials",
        "FodderGrasses",
        "Rye",
        "Barley",
        "Millet",
        "Sunflower",
        "GroundnutsPeanuts",
        "RapeseedCanola",
        "Potatoes",
        "Sugarbeet",
        "Cassava",
    ]

    # jianyong01b
    this = "jianyong01b"
    crops = ["Wheat", "Maize", "Sorghum", "Rice", "Soybean", "FabaBean"]
    out_keys[this] = [["Wheat"], ["Maize"], ["Sorghum"], ["Rice"], ["Soybean"], ["Pulses"]]
    out_list[this] = crops
    out_ignores[this] = [
        "Sugar cane",
        "Oil palm",
        "Citrus",
        "Date palm",
        "Grapes / Vine",
        "Cotton",
        "Cocoa",
        "Coffee",
        "Others annual",
        "Others perennial",
        "Fodder grasses",
        "Rye",
        "Barley",
        "Millet",
        "Sunflower",
        "Groundnuts / Peanuts",
        "Rape seed / Canola",
        "Potatoes",
        "Sugar beet",
        "Cassava",
    ]

    # ani01
    this = "ani01"
    crops = ["TeWW", "TeCo", "TrRi", "TeSW"]
    out_keys[this] = [
        ["Wheat", "Rye", "Barley"],
        ["Maize", "Millet", "Sorghum", "Sugar cane"],
        ["Rice"],
        [
            "Soybean",
            "Pulses",
            "Oil palm",
            "Citrus",
            "Date palm",
            "Grapes / Vine",
            "Cotton",
            "Cocoa",
            "Coffee",
            "Others annual",
            "Others perennial",
            "Fodder grasses",
            "Sunflower",
            "Groundnuts / Peanuts",
            "Rape seed / Canola",
            "Potatoes",
            "Sugar beet",
            "Cassava",
        ],
    ]
    out_list[this] = crops
    out_ignores[this] = []

    # 20200928
    this = "20200928"
    crops = ["whe", "mai", "ric", "soy", "bea", "sor", "mil"]
    out_keys[this] = [["Wheat"], ["Maize"], ["Rice"], ["Soybeans"], ["Pulses"], ["Sorghum"], ["Millet"]]
    out_list[this] = crops
    out_ignores[this] = [
        "Sugarcane",
        "Citrus",
        "Datepalm",
        "GrapesVine",
        "Cotton",
        "Cocoa",
        "Coffee",
        "OtherAnnuals",
        "OtherPerennials",
        "FodderGrasses",
        "Sugarbeet",
        "Barley",
        "Rye",
        "Sunflower",
        "GroundnutsPeanuts",
        "RapeseedCanola",
        "Oilpalm",
        "Potatoes",
        "Cassava",
    ]

    # ignore types grouping
    all_ignores = [out_ignores[n] for n in out_names]
    out_ignore_types = [0] * len(all_ignores)
    for i in range(len(all_ignores)):
        if out_ignore_types[i] == 0:
            out_ignore_types[i] = i + 1
            this_ign = all_ignores[i]
            for j in range(i + 1, len(all_ignores)):
                if this_ign == all_ignores[j]:
                    out_ignore_types[j] = i + 1

    list_crops_out = out_list[this_ver]
    in2out_key = out_keys[this_ver]
    list_ignore = out_ignores[this_ver]

    return (
        list_crops_out,
        in2out_key,
        list_ignore,
        out_names,
        [out_keys[n] for n in out_names],
        all_ignores,
        out_ignore_types,
    )


def process_mirca(
    mirca_path: str,
    gridlist,
    this_ver: str,
    force_all_rainfed: bool,
    plum_setaside_frac: float,
):
    """Read MIRCA data and map to output crops (MATLAB remap.m)."""
    print("Importing crop fractions...")
    print(f"file_cropmirca: {mirca_path}")
    croparea_in = read_table_then2map(mirca_path, verboseIfNoMat=True)
    tmp = yxz_to_xz(croparea_in["maps_YXv"], (len(gridlist["list2map"]), croparea_in["maps_YXv"].shape[2]), gridlist["list2map"])
    croparea_in.pop("maps_YXv")
    croparea_in["garr_xv"] = tmp

    list_crops_frac_in = list(croparea_in["varNames"])
    list_crops_combined_in = sorted({c.replace("_IR", "").replace("_RF", "") for c in list_crops_frac_in})

    if force_all_rainfed:
        print("Warning: GETTING RID OF IRRIGATED")
        for crop in list_crops_combined_in:
            ir = croparea_in["varNames"].index(f"{crop}_IR")
            rf = croparea_in["varNames"].index(f"{crop}_RF")
            croparea_in["garr_xv"][:, rf] += croparea_in["garr_xv"][:, ir]
            croparea_in["garr_xv"][:, ir] = 0 * croparea_in["garr_xv"][:, ir]

    for crop in list_crops_combined_in:
        is_this = [crop in name for name in croparea_in["varNames"]]
        if sum(is_this) != 2:
            raise RuntimeError("length(find(isThisCrop))~=2")

    list_crops_combined_out, in2out_key_combined, list_ignore_frac, all_names, _, _, all_ignore_types = get_remapv2_keys(this_ver)
    list_crops_out = list_crops_combined_out + [f"{c}i" for c in list_crops_combined_out]

    if len(list_crops_combined_out) != len(set(list_crops_combined_out)):
        raise RuntimeError("length(list_cropsCombined_out) ~= length(unique(list_cropsCombined_out))")

    # Check that every crop in croparea_in is included or ignored
    for crop in list_crops_combined_in:
        found = 1 if crop in list_ignore_frac else 0
        for lst in in2out_key_combined:
            found += len([x for x in lst if x == crop])
        if found == 0:
            raise RuntimeError(f"{crop} not found in list_ignore_frac or in2out_keyCombined_frac!")
        if found > 1:
            raise RuntimeError(f"{crop} found {found} times in list_ignore_frac and in2out_keyCombined_frac!")

    if any(len(x) == 0 for x in in2out_key_combined):
        raise RuntimeError("At least one member of in2out_keyCombined_frac is empty!")

    print("Crop fraction mappings (out: in):")
    for c, comp in enumerate(in2out_key_combined):
        print(f"    {list_crops_combined_out[c]}: {' '.join(comp)}")

    in2out_key_frac: list[list[str]] = []
    for c in range(len(list_crops_combined_out)):
        in2out_key_frac.append([f"{x}_RF" for x in in2out_key_combined[c]])
    for c in range(len(list_crops_combined_out)):
        in2out_key_frac.append([f"{x}_IR" for x in in2out_key_combined[c]])

    n_cells = len(gridlist["list2map"])
    croparea_mid = {"garr_xv": np.full((n_cells, len(list_crops_out)), np.nan), "varNames": list_crops_out}
    for c, this_crop in enumerate(list_crops_out):
        this_row = in2out_key_frac[c]
        ia = [croparea_in["varNames"].index(x) for x in this_row]
        if len(ia) != len(this_row):
            raise RuntimeError("length(IA)~=length(thisRow)")
        croparea_mid["garr_xv"][:, c] = np.sum(croparea_in["garr_xv"][:, ia], axis=1)

    # Move ignored area and setAside area into ExtraCrop
    ignore_mask = [any(ign in name for ign in list_ignore_frac) for name in croparea_in["varNames"]]
    ignore_area_x = np.sum(croparea_in["garr_xv"][:, np.array(ignore_mask)], axis=1)
    setaside_area_x = np.sum(croparea_mid["garr_xv"] * plum_setaside_frac, axis=1)
    croparea_mid["garr_xv"] = croparea_mid["garr_xv"] * (1 - plum_setaside_frac)
    extra_x = ignore_area_x + setaside_area_x
    croparea_mid["garr_xv"] = np.column_stack([croparea_mid["garr_xv"], extra_x])
    list_crops_combined_out = list_crops_combined_out + ["ExtraCrop"]
    list_crops_out = list_crops_out + ["ExtraCrop"]

    mid_cropfrac = {
        "garr_xv": croparea_mid["garr_xv"]
        / np.tile(np.sum(croparea_mid["garr_xv"], axis=1).reshape(-1, 1), (1, len(list_crops_out))),
        "varNames": list_crops_out,
        "list2map": gridlist["list2map"],
        "lonlats": gridlist["lonlats"],
    }

    print("Done.")
    return (
        mid_cropfrac,
        croparea_in,
        list_crops_out,
        list_crops_combined_out,
        list_ignore_frac,
        in2out_key_frac,
        all_names,
        all_ignore_types,
    )
