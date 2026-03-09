"""Standard entrypoint for PLUMharmFigs diagnostic figures (all scenarios).

Equivalent to running PLUMharm_options.m + PLUMharmFigs_options.m then PLUMharmFigs.m
in MATLAB.
"""
from __future__ import annotations

from pathlib import Path

from landsymm.config import get_project_root, get_remap_output_dir
from .plumharm_figs import run_plumharm_figs
from .plumharm_figs_options import PlumHarmFigsConfig

SCENARIOS = [
    "SSP1_RCP26",
    "SSP2_RCP45",
    "SSP3_RCP70",
    "SSP4_RCP60",
    "SSP5_RCP85",
]


def main(scenarios: list[str] | None = None, year1: int = 2021, yearN: int = 2100) -> None:
    repo_root = get_project_root()
    data_dir = repo_root / "data" / "PLUMv2_LU_default_output"
    baseline_dir = get_remap_output_dir() / "remaps_v10_old_62892_gL"

    if scenarios is None:
        scenarios = SCENARIOS

    for ssp in scenarios:
        print(f"\n{'='*60}")
        print(f"  Running PLUMharmFigs for {ssp}")
        print(f"{'='*60}")

        harm_dir = str(data_dir / ssp / "s1.HILDA+_remap_v10_old_62892_gL.harm.allow_unveg_py")
        figs_dir = str(data_dir / ssp / "s1.HILDA+_remap_v10_old_62892_gL.harm.allow_unveg_py_harms_figs")

        cfg = PlumHarmFigsConfig(
            this_dir=str(repo_root),
            geodata_dir=str(repo_root / "data" / "geodata_py"),
            plum_dirs=[str(data_dir / ssp / "s1")],
            harm_dirs=[harm_dir],
            remap_lu_file=str(baseline_dir / "LU.remapv10_old_62892_gL.txt"),
            remap_cropf_file=str(baseline_dir / "cropfracs.remapv10_old_62892_gL.txt"),
            remap_nfert_file=str(baseline_dir / "nfert.remapv10_old_62892_gL.txt"),
            base_year=2020,
            year1=year1,
            yearN=yearN,
            combine_crops=False,
            fruitveg_sugar_2oil=False,
            allow_unveg=True,
            norm2extra=0.177,
            timeseries_legend_loc="best",
            harms_figs_dir=figs_dir,
            runlist_legend=[ssp],
            save_geotiffs=False,
            this_ver="",
            three_years=[year1, (year1 + yearN) // 2, yearN],
            year_list_baseline_to_plot=None,
        )

        run_plumharm_figs(cfg)

    print("\nAll scenarios complete.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run PLUMharmFigs for all scenarios.")
    parser.add_argument(
        "--scenarios", nargs="+", default=None,
        help=f"Scenarios to run (default: all). Options: {SCENARIOS}",
    )
    parser.add_argument("--year1", type=int, default=2021, help="First year (default: 2021)")
    parser.add_argument("--yearN", type=int, default=2100, help="Last year (default: 2100)")
    args = parser.parse_args()

    main(scenarios=args.scenarios, year1=args.year1, yearN=args.yearN)
