"""Run the full PLUM harmonization pipeline (all scenarios).

Executes three stages in sequence:
  1. Reformat raw PLUM outputs (reformat_plum_gridded)
  2. PLUMharm harmonization (run_plumharm)
  3. PLUMharm2LPJG conversion to LPJ-GUESS inputs (run_plumharm2lpjg)

Equivalent to running reformat_gridded_updated.R, PLUMharm.m,
then PLUMharm2LPJG.m in the MATLAB/R workflow.
"""
from __future__ import annotations

import time

from landsymm.config import discover_scenarios, get_member, get_plum_output_dir
from .reformat_plum_gridded import reformat_all
from .run_plumharm import main as run_harm
from .run_plumharm2lpjg import main as run_2lpjg


def main(
    scenarios: list[str] | None = None,
    year1: int = 2021,
    yearN: int = 2100,
    parent_dir: str | None = None,
    skip_reformat: bool = False,
) -> None:
    t0 = time.perf_counter()

    if parent_dir is None:
        parent_dir = str(get_plum_output_dir())

    shown = scenarios or discover_scenarios(parent_dir, get_member()) or "(none found)"
    print("=" * 60)
    print("  PLUM Harmonization Pipeline")
    print(f"  Scenarios: {shown}")
    print(f"  Years: {year1}-{yearN}")
    print(f"  Parent dir: {parent_dir}")
    print("=" * 60)

    if not skip_reformat:
        print("\n>>> Step 1: Reformat raw PLUM outputs")
        reformat_all(parent_dir, scenarios=scenarios)
    else:
        print("\n>>> Step 1: Reformat (SKIPPED)")

    print("\n>>> Step 2: PLUMharm (harmonization)")
    run_harm(scenarios=scenarios, year1=year1, yearN=yearN, parent_dir=parent_dir)

    print("\n>>> Step 3: PLUMharm2LPJG (conversion to LPJ-GUESS inputs)")
    run_2lpjg(scenarios=scenarios, year1=year1, yearN=yearN, parent_dir=parent_dir)

    elapsed = time.perf_counter() - t0
    minutes = int(elapsed // 60)
    seconds = int(elapsed % 60)
    print(f"\n{'='*60}")
    print(f"  Pipeline complete ({minutes}m {seconds}s)")
    print(f"{'='*60}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Run full PLUM harmonization pipeline (reformat + PLUMharm + PLUMharm2LPJG)."
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
    parser.add_argument(
        "--skip-reformat", action="store_true",
        help="Skip the reformatting step (if already done)",
    )
    args = parser.parse_args()

    main(
        scenarios=args.scenarios,
        year1=args.year1,
        yearN=args.yearN,
        parent_dir=args.parent_dir,
        skip_reformat=args.skip_reformat,
    )
