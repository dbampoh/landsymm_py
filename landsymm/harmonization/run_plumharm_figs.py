"""Standard entrypoint for PLUMharmFigs diagnostic figures (all scenarios).

Equivalent to running PLUMharm_options.m + PLUMharmFigs_options.m then PLUMharmFigs.m
in MATLAB. Scenarios default to all found under the PLUM parent dir; data locations
and naming are configurable via the ``LANDSYMM_*`` environment variables documented
in ``landsymm.config``.
"""
from __future__ import annotations

from pathlib import Path

from landsymm.config import (
    discover_scenarios,
    get_geodata_dir,
    get_member,
    get_plum_output_dir,
    get_project_root,
    get_remap_baseline_files,
    harm_dirname,
)
from .plumharm_figs import run_plumharm_figs
from .plumharm_figs_options import PlumHarmFigsConfig


def main(
    scenarios: list[str] | None = None,
    year1: int = 2021,
    yearN: int = 2100,
    parent_dir: str | None = None,
) -> None:
    repo_root = get_project_root()
    data_dir = Path(parent_dir) if parent_dir else get_plum_output_dir()
    member = get_member()
    baseline = get_remap_baseline_files()
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
        print(f"  Running PLUMharmFigs for {ssp}")
        print(f"{'='*60}")

        harm_dir = str(data_dir / ssp / harm_sub)
        figs_dir = str(data_dir / ssp / f"{harm_sub}_harms_figs")

        cfg = PlumHarmFigsConfig(
            this_dir=str(repo_root),
            geodata_dir=str(get_geodata_dir()),
            plum_dirs=[str(data_dir / ssp / member)],
            harm_dirs=[harm_dir],
            remap_lu_file=str(baseline["lu"]),
            remap_cropf_file=str(baseline["cropf"]),
            remap_nfert_file=str(baseline["nfert"]),
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

    parser = argparse.ArgumentParser(
        description="Run PLUMharmFigs (default: all scenarios found under the PLUM parent dir)."
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
