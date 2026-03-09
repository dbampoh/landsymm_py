"""Run the full wetland-into-HILDA pipeline.

Step 1: Aggregate GLWD3 30-arcsec raster to half-degree wetland fractions.
Step 2: Insert GLWD3 wetland/peatland into HILDA+ remap LU data (Approach H).
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

from landsymm.config import get_project_root
from .glwd3_to_halfdeg import aggregate_glwd3_to_halfdeg
from .wetland_into_hilda import insert_wetland_approach_h


def main(
    glwd3_path: str | None = None,
    lu_path: str | None = None,
    output_path: str | None = None,
    glwd3_output_dir: str | None = None,
    wetland_product: str = "wforests",
    skip_aggregation: bool = False,
    verbose: bool = True,
) -> None:
    t0 = time.perf_counter()

    repo_root = get_project_root()
    glwd3_dir = glwd3_output_dir or str(repo_root / "data" / "geodata_py" / "glwd3")

    if glwd3_path is None:
        glwd3_path = str(Path(glwd3_dir) / "glwd_3.tif")
    if lu_path is None:
        from landsymm.config import get_remap_output_dir
        lu_path = str(
            get_remap_output_dir()
            / "remaps_v10_old_62892_gL" / "LU.remapv10_old_62892_gL.txt"
        )

    wetland_nc = str(Path(glwd3_dir) / "peatland_halfdeg.nc")

    if output_path is None:
        import os
        lu_dir = os.path.dirname(lu_path)
        lu_base = os.path.basename(lu_path)
        name, ext = os.path.splitext(lu_base)
        output_path = os.path.join(lu_dir, f"{name}_peatland{ext}")

    print("=" * 60)
    print("  Wetland-into-HILDA Pipeline")
    print(f"  GLWD3 source: {glwd3_path}")
    print(f"  LU input: {lu_path}")
    print(f"  LU output: {output_path}")
    print(f"  Wetland product: {wetland_product}")
    print("=" * 60)

    # Step 1
    if not skip_aggregation:
        print("\n>>> Step 1: Aggregate GLWD3 to half-degree wetland fractions")
        aggregate_glwd3_to_halfdeg(glwd3_path, glwd3_dir, verbose=verbose)
    else:
        print("\n>>> Step 1: Aggregation (SKIPPED)")

    # Step 2
    print("\n>>> Step 2: Insert wetland into HILDA+ remap LU (Approach H)")
    insert_wetland_approach_h(
        lu_path, wetland_nc, output_path,
        wetland_product=wetland_product, verbose=verbose,
    )

    elapsed = time.perf_counter() - t0
    minutes = int(elapsed // 60)
    seconds = int(elapsed % 60)
    print(f"\n{'=' * 60}")
    print(f"  Pipeline complete ({minutes}m {seconds}s)")
    print(f"  Output: {output_path}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run full wetland-into-HILDA pipeline (GLWD3 aggregation + Approach H insertion)."
    )
    parser.add_argument(
        "--glwd3-path", default=None,
        help="Path to GLWD3 GeoTIFF (default: data/geodata_py/glwd3/glwd_3.tif)",
    )
    parser.add_argument(
        "--lu-path", default=None,
        help="Path to HILDA+ remap LU file",
    )
    parser.add_argument(
        "--output-path", default=None,
        help="Output LU file path (default: alongside input with _peatland tag)",
    )
    parser.add_argument(
        "--glwd3-output-dir", default=None,
        help="Directory for GLWD3 intermediate outputs (default: data/geodata_py/glwd3/)",
    )
    parser.add_argument(
        "--wetland-product",
        choices=["wforests", "noforests"],
        default="wforests",
        help=(
            "Which wetland product to use: 'wforests' includes forest wetland "
            "types (classes 5,6); 'noforests' excludes them (default: wforests)"
        ),
    )
    parser.add_argument(
        "--skip-aggregation", action="store_true",
        help="Skip GLWD3 aggregation (if peatland_halfdeg.nc already exists)",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    main(
        glwd3_path=args.glwd3_path,
        lu_path=args.lu_path,
        output_path=args.output_path,
        glwd3_output_dir=args.glwd3_output_dir,
        wetland_product=args.wetland_product,
        skip_aggregation=args.skip_aggregation,
        verbose=not args.quiet,
    )
