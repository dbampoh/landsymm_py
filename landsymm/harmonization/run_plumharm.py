"""Standard entrypoint for PLUM harmonization (all scenarios).

Equivalent to running PLUMharm_options.m then PLUMharm.m in MATLAB.
"""
from __future__ import annotations

from pathlib import Path

from landsymm.config import get_project_root, get_remap_output_dir
from .plumharm import run_plumharm
from .plumharm_options import PlumHarmConfig

SCENARIOS = [
    "SSP1_RCP26",
    "SSP2_RCP45",
    "SSP3_RCP70",
    "SSP4_RCP60",
    "SSP5_RCP85",
]


def main(
    scenarios: list[str] | None = None,
    year1: int = 2021,
    yearN: int = 2100,
    parent_dir: str | None = None,
) -> None:
    repo_root = get_project_root()
    data_dir = Path(parent_dir) if parent_dir else repo_root / "data" / "PLUMv2_LU_default_output"
    baseline_dir = get_remap_output_dir() / "remaps_v10_old_62892_gL"

    if scenarios is None:
        scenarios = SCENARIOS

    for ssp in scenarios:
        print(f"\n{'='*60}")
        print(f"  Running PLUMharm for {ssp}")
        print(f"{'='*60}")

        cfg = PlumHarmConfig(
            this_dir=str(repo_root),
            geodata_dir=str(repo_root / "data" / "geodata_py"),
            plum_dirs=[str(data_dir / ssp / "s1")],
            harm_dirs=[
                str(data_dir / ssp / f"s1.HILDA+_remap_v10_old_62892_gL.harm.allow_unveg_py")
            ],
            remap_lu_file=str(baseline_dir / "LU.remapv10_old_62892_gL.txt"),
            remap_cropf_file=str(baseline_dir / "cropfracs.remapv10_old_62892_gL.txt"),
            remap_nfert_file=str(baseline_dir / "nfert.remapv10_old_62892_gL.txt"),
            base_year=2020,
            year1=year1,
            yearN=yearN,
            allow_unveg=True,
            conserv_tol_pct=0.2,
            conserv_tol_area=1e3,
            norm2extra=0.177,
            use_latest_plum_mgmt=True,
            inpaint_method=4,
            fix_tiny_negs_tol_m2=1.0,
            save_halfdeg_mat=True,
            save_2deg_mat=True,
            save_halfdeg_txt=False,
            save_2deg_txt=False,
            out_prec=6,
            out_width=1,
            delimiter=" ",
            overwrite=True,
            fancy=False,
            debug_areas=False,
            debug_nfert=False,
            debug_irrig=False,
            debugIJ_2deg=None,
            dbCrop="",
            verbose=True,
            combine_crops=False,
            fruitveg_sugar_2oil=False,
        )

        run_plumharm(cfg)

    print("\nAll scenarios complete.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run PLUMharm for all scenarios.")
    parser.add_argument(
        "--scenarios", nargs="+", default=None,
        help=f"Scenarios to run (default: all). Options: {SCENARIOS}",
    )
    parser.add_argument("--year1", type=int, default=2021, help="First year (default: 2021)")
    parser.add_argument("--yearN", type=int, default=2100, help="Last year (default: 2100)")
    parser.add_argument(
        "--parent-dir", default=None,
        help="Parent directory containing SSP scenario dirs (default: data/PLUMv2_LU_default_output)",
    )
    args = parser.parse_args()

    main(scenarios=args.scenarios, year1=args.year1, yearN=args.yearN, parent_dir=args.parent_dir)
