from landsymm.config import get_geodata_dir, get_remap_output_dir
from landsymm.remapping.remap import run_remap
from landsymm.remapping.remap_options import RemapConfig


def main() -> None:
    geodata = get_geodata_dir()
    gridlists = geodata / "gridlists"
    soil = geodata / "soil"

    cfg = RemapConfig(
        geodata_dir=str(geodata),
        year_list_out=list(range(1901, 2021)),
        lu_source="HILDA+",
        out_dir_top=str(get_remap_output_dir()),
        fill_unveg=1e-6,
        interp_test_y1=1995,
        interp_test_yN=2005,
        inpaint_method=4,
        force_all_rainfed=False,
        plum_setaside_frac=0.088,
        remap_ver="10_old_62892_gL",
        this_ver="WithFruitVeg_sepOil_combSugar",
        file_gridlist_out=str(gridlists / "gridlist_62892.runAEclimOK.txt"),
        file_gridlist_climate=str(
            gridlists / "remaps_v10_g2p_isimipclimMask"
            / "gridlist.remapv10_g2p_isimipclimMask.txt"
        ),
        files_soil=(
            str(soil / "HWSD_soil_data_on_cropland_v2.2.remapv10_old_62892_gL.txt"),
            str(soil / "HWSD_soil_data_on_cropland_v2.3.remapv10_old_62892_gL.txt"),
            str(soil / "soilmap_center_interpolated.dat"),
        ),
        out_prec=6,
        out_width=1,
        delimiter=" ",
        overwrite=True,
        fancy=False,
        do_gzip=False,
        log_file=None,
    )

    run_remap(cfg)


if __name__ == "__main__":
    main()
