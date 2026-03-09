#!/usr/bin/env python3
"""
Inspect upscaled netfrac output and optionally validate forestfrac output.
"""

import argparse
import random
from pathlib import Path

import numpy as np

DEFAULT_OUTPUT_DIR = str(Path(__file__).resolve().parents[3] / 'data' / 'geodata_py' / 'HILDA+' / 'data' / 'output')
DEFAULT_REPORT_NAME = "inspect_netfrac_forestfrac.txt"


def parse_args():
    parser = argparse.ArgumentParser(description="Inspect netfrac output")
    parser.add_argument("--input", required=True, help="Path to netfrac output file")
    parser.add_argument("--forestfrac", default=None,
                        help="Optional forestfrac output file to validate")
    parser.add_argument("--gridlist", default=None,
                        help="Optional gridlist file to compare gridcell counts")
    parser.add_argument(
        "--output",
        default=None,
        help=f"Path to write report (default: {DEFAULT_OUTPUT_DIR}/{DEFAULT_REPORT_NAME})",
    )
    parser.add_argument("--sample-rows", type=int, default=5,
                        help="Number of sample rows to include (default: 5)")
    parser.add_argument("--sum-check-rows", type=int, default=500,
                        help="Rows to use for sum-to-1 check (default: 500)")
    parser.add_argument("--random-samples", type=int, default=0,
                        help="Number of random rows to include (default: 0)")
    parser.add_argument("--random-seed", type=int, default=None,
                        help="Random seed for reproducible sampling")
    parser.add_argument("--forest-sum-tol", type=float, default=1e-6,
                        help="Tolerance for forestfrac sum checks (default: 1e-6)")
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output) if args.output else Path(DEFAULT_OUTPUT_DIR) / DEFAULT_REPORT_NAME

    with open(input_path, "r") as f:
        header = f.readline().rstrip("\n")
    net_header_fields = header.split()
    forest_idx = net_header_fields.index("FOREST") if "FOREST" in net_header_fields else None

    report = []
    report.append(f"File: {input_path}")
    report.append(f"Header: {header}")

    sample = np.loadtxt(input_path, skiprows=1, max_rows=args.sample_rows)
    report.append("Sample rows:")
    report.append(np.array2string(sample, max_line_width=120))

    subset = np.loadtxt(input_path, skiprows=1, max_rows=args.sum_check_rows)
    fractions = subset[:, 3:]
    row_sums = fractions.sum(axis=1)
    report.append(f"Sum check (first 5): {row_sums[:5].tolist()}")
    report.append(f"Sum min/max (first {args.sum_check_rows}): {row_sums.min()} / {row_sums.max()}")

    unique_gridcells = set()
    unique_years = set()
    with open(input_path, "r") as f:
        f.readline()
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            unique_gridcells.add((parts[0], parts[1]))
            unique_years.add(parts[2])
    report.append(f"Unique gridcells: {len(unique_gridcells)}")
    if unique_years:
        report.append(f"Unique years: {len(unique_years)} (range: {min(unique_years)}-{max(unique_years)})")
    if args.gridlist:
        gridlist_cells = set()
        with open(args.gridlist, "r") as f:
            for line in f:
                if not line.strip():
                    continue
                if line.lstrip().startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                gridlist_cells.add((parts[0], parts[1]))
        report.append(f"Gridlist gridcells: {len(gridlist_cells)}")

    if args.random_samples and args.random_samples > 0 and not args.forestfrac:
        data = np.loadtxt(input_path, skiprows=1)
        nrows = data.shape[0]
        if args.random_seed is not None:
            np.random.seed(args.random_seed)
        idx = np.random.choice(nrows, size=min(args.random_samples, nrows), replace=False)
        report.append(f"Random sample rows (count={len(idx)}):")
        report.append(np.array2string(data[idx], max_line_width=120))

    if args.forestfrac:
        forest_path = Path(args.forestfrac)
        report.append("")
        report.append(f"Forestfrac file: {forest_path}")

        rng = random.Random(args.random_seed)
        sample_size = max(0, int(args.random_samples))
        sample_rows = []
        total_rows = 0
        mismatch_count = 0
        mismatch_example = None

        sum_min = float("inf")
        sum_max = float("-inf")
        sum_min_nonzero = float("inf")
        sum_max_nonzero = float("-inf")
        count_nonzero_forest = 0
        count_zero_forest = 0
        count_nonzero_bad = 0
        count_zero_bad = 0
        forest_abs_diff_max = 0.0

        with open(input_path, "r") as nf, open(forest_path, "r") as ff:
            nf_header = nf.readline().strip()
            ff_header = ff.readline().strip()
            report.append(f"Forestfrac header: {ff_header}")

            for nline, fline in zip(nf, ff):
                nline = nline.strip()
                fline = fline.strip()
                if not nline or not fline:
                    continue

                total_rows += 1
                nvals = nline.split()
                fvals = fline.split()

                if nvals[:3] != fvals[:3]:
                    mismatch_count += 1
                    if mismatch_example is None:
                        mismatch_example = (nvals[:3], fvals[:3])

                fdata = list(map(float, fvals[3:]))
                fsum = sum(fdata)
                sum_min = min(sum_min, fsum)
                sum_max = max(sum_max, fsum)

                forest_val = float(nvals[forest_idx]) if forest_idx is not None else None
                if forest_val is not None:
                    forest_abs_diff = abs((fsum * forest_val) - forest_val)
                    forest_abs_diff_max = max(forest_abs_diff_max, forest_abs_diff)
                    if forest_val > 0:
                        count_nonzero_forest += 1
                        sum_min_nonzero = min(sum_min_nonzero, fsum)
                        sum_max_nonzero = max(sum_max_nonzero, fsum)
                        if abs(fsum - 1.0) > args.forest_sum_tol:
                            count_nonzero_bad += 1
                    else:
                        count_zero_forest += 1
                        if abs(fsum) > args.forest_sum_tol:
                            count_zero_bad += 1

                if sample_size > 0:
                    if len(sample_rows) < sample_size:
                        sample_rows.append((nvals[:3], forest_val, fsum))
                    else:
                        j = rng.randint(1, total_rows)
                        if j <= sample_size:
                            sample_rows[j - 1] = (nvals[:3], forest_val, fsum)

        report.append(f"Total rows checked: {total_rows}")
        report.append(f"Row alignment mismatches: {mismatch_count}")
        if mismatch_example:
            report.append(f"First mismatch (netfrac vs forestfrac): {mismatch_example}")
        report.append(f"Forestfrac sum min/max (all): {sum_min} / {sum_max}")
        if forest_idx is not None and count_nonzero_forest:
            report.append(f"Forestfrac sum min/max (forest>0): {sum_min_nonzero} / {sum_max_nonzero}")
            report.append(f"Rows with forest>0: {count_nonzero_forest}, bad sums: {count_nonzero_bad}")
            report.append(f"Rows with forest=0: {count_zero_forest}, nonzero sums: {count_zero_bad}")
            report.append(f"Max |FOREST*(sum forestfrac) - FOREST|: {forest_abs_diff_max}")
        elif forest_idx is None:
            report.append("Note: FOREST column not found in netfrac header; skipped forest>0 checks.")

        if sample_rows:
            report.append(f"Random sample rows (count={len(sample_rows)}) [lon lat year, FOREST, forestfrac sum]:")
            report.append(np.array2string(np.array(sample_rows, dtype=object), max_line_width=120))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(report) + "\n")
    print(f"Wrote report to: {output_path}")


if __name__ == "__main__":
    main()
