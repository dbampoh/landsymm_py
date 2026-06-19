"""Standard entrypoint for PLUMharm2LPJG conversion (all scenarios).

Equivalent to running PLUMharm_options.m + PLUMharm2LPJG_options.m then PLUMharm2LPJG.m
in MATLAB.

Scenarios default to every scenario directory found under the PLUM parent dir.
Data locations and naming are configurable via the ``LANDSYMM_*`` environment
variables documented in ``landsymm.config``.
"""
from __future__ import annotations

from pathlib import Path

from landsymm.config import (
    discover_scenarios,
    get_member,
    get_plum_output_dir,
    harm_dirname,
)
from .plumharm2lpjg import run_plumharm2lpjg
from .plumharm2lpjg_options import PlumHarm2LPJGConfig


def main(
    scenarios: list[str] | None = None,
    year1: int = 2021,
    yearN: int = 2100,
    parent_dir: str | None = None,
) -> None:
    data_dir = Path(parent_dir) if parent_dir else get_plum_output_dir()
    member = get_member()
    harm_sub = harm_dirname(member=member, allow_unveg=True)

    if scenarios is None:
        scenarios = discover_scenarios(data_dir, member)
    if not scenarios:
        raise RuntimeError(
            f"No scenarios found under {data_dir} "
            f"(expected subdirectories each containing a '{member}' dir)."
        )

    for ssp in scenarios:
        print(f"\n{'='*60}")
        print(f"  Running PLUMharm2LPJG for {ssp}")
        print(f"{'='*60}")

        cfg = PlumHarm2LPJGConfig(
            this_dir=str(data_dir),
            plum_dirs=[f"{ssp}/{member}"],
            harm_dirs=[f"{ssp}/{harm_sub}"],
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

    parser = argparse.ArgumentParser(
        description="Run PLUMharm2LPJG (default: all scenarios found under the PLUM parent dir)."
    )
    parser.add_argument(
        "--scenarios", nargs="+", default=None,
        help="Scenarios to run (default: all discovered under the parent dir).",
    )
    parser.add_argument("--year1", type=int, default=2021, help="First year (default: 2021)")
    parser.add_argument("--yearN", type=int, default=2100, help="Last year (default: 2100)")
    parser.add_argument(
        "--parent-dir", default=None,
        help="PLUM scenario parent dir (default: config.get_plum_output_dir()).",
    )
    args = parser.parse_args()

    main(scenarios=args.scenarios, year1=args.year1, yearN=args.yearN, parent_dir=args.parent_dir)
