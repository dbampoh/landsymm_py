"""Standard entrypoint for PLUMharm2LPJG conversion (all scenarios).

Equivalent to running PLUMharm_options.m + PLUMharm2LPJG_options.m then PLUMharm2LPJG.m
in MATLAB.
"""
from __future__ import annotations

from pathlib import Path

from landsymm.config import get_project_root
from .plumharm2lpjg import run_plumharm2lpjg
from .plumharm2lpjg_options import PlumHarm2LPJGConfig

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

    if scenarios is None:
        scenarios = SCENARIOS

    for ssp in scenarios:
        print(f"\n{'='*60}")
        print(f"  Running PLUMharm2LPJG for {ssp}")
        print(f"{'='*60}")

        cfg = PlumHarm2LPJGConfig(
            this_dir=str(data_dir),
            plum_dirs=[f"{ssp}/s1"],
            harm_dirs=[
                f"{ssp}/s1.HILDA+_remap_v10_old_62892_gL.harm.allow_unveg_py"
            ],
            base_year=2020,
            year1=year1,
            yearN=yearN,
            combine_crops=False,
            fruitveg_sugar_2oil=False,
            out_prec=6,
            out_width=1,
            delimiter=" ",
            overwrite=True,
            fancy=False,
            do_gzip=False,
            donation_order=("PASTURE", "NATURAL", "BARREN"),
            save_every_pct=1,
            someofall=True,
            forLPJG_dirs=None,
            verbose_write=False,
            yStep=1,
            y1_pre=None,
        )

        run_plumharm2lpjg(cfg)

    print("\nAll scenarios complete.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run PLUMharm2LPJG for all scenarios.")
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
