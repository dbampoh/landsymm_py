"""Standard entrypoint for PLUM harmonization (all scenarios).

Equivalent to running PLUMharm_options.m then PLUMharm.m in MATLAB.

Scenarios default to every scenario directory found under the PLUM parent dir
(see ``config.discover_scenarios``). Data locations and naming are configurable
via the ``LANDSYMM_*`` environment variables documented in ``landsymm.config``
(no source edits needed for a new scenario set or data location).
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
from .plumharm import run_plumharm
from .plumharm_options import PlumHarmConfig


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
        print(f"  Running PLUMharm for {ssp}")
        print(f"{'='*60}")

        cfg = PlumHarmConfig(
            this_dir=str(repo_root),
            geodata_dir=str(get_geodata_dir()),
            plum_dirs=[str(data_dir / ssp / member)],
            harm_dirs=[str(data_dir / ssp / harm_sub)],
            remap_lu_file=str(baseline["lu"]),
            remap_cropf_file=str(baseline["cropf"]),
            remap_nfert_file=str(baseline["nfert"]),
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

    parser = argparse.ArgumentParser(
        description="Run PLUMharm (default: all scenarios found under the PLUM parent dir)."
    )
    parser.add_argument(
        "--scenarios", nargs="+", default=None,
        help="Scenarios to run (default: all discovered under the parent dir).",
    )
    parser.add_argument("--year1", type=int, default=2021, help="First year (default: 2021)")
    parser.add_argument("--yearN", type=int, default=2100, help="Last year (default: 2100)")
    parser.add_argument(
        "--parent-dir", default=None,
        help="PLUM scenario parent dir (default: config.get_plum_output_dir(), "
             "set via LANDSYMM_DATA_DIR / LANDSYMM_PLUM_DIRNAME).",
    )
    args = parser.parse_args()

    main(scenarios=args.scenarios, year1=args.year1, yearN=args.yearN, parent_dir=args.parent_dir)
