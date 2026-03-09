#!/usr/bin/env python3
"""
Inspect a smoothed HILDA+ NetCDF file and write a brief report.
Works for full or mini outputs.
"""

import argparse
from pathlib import Path

import netCDF4 as nc
import numpy as np

DEFAULT_OUTPUT_DIR = str(Path(__file__).resolve().parents[3] / 'data' / 'geodata_py' / 'HILDA+' / 'data' / 'output')


def parse_args():
    parser = argparse.ArgumentParser(description="Inspect smoothed HILDA+ NetCDF")
    parser.add_argument("--input", required=True, help="Path to smoothed NetCDF file")
    parser.add_argument(
        "--output",
        default=None,
        help=f"Path to write report (default: {DEFAULT_OUTPUT_DIR}/inspect_smoothed.txt)",
    )
    parser.add_argument("--sample-lat", type=int, default=None,
                        help="Starting latitude index for sample window")
    parser.add_argument("--sample-lon", type=int, default=None,
                        help="Starting longitude index for sample window")
    parser.add_argument("--sample-lat-min", type=int, default=None,
                        help="Minimum latitude index for sample window")
    parser.add_argument("--sample-lat-max", type=int, default=None,
                        help="Maximum latitude index for sample window (exclusive)")
    parser.add_argument("--sample-lon-min", type=int, default=None,
                        help="Minimum longitude index for sample window")
    parser.add_argument("--sample-lon-max", type=int, default=None,
                        help="Maximum longitude index for sample window (exclusive)")
    parser.add_argument("--sample-lat-min-deg", type=float, default=None,
                        help="Minimum latitude in degrees for sample window")
    parser.add_argument("--sample-lat-max-deg", type=float, default=None,
                        help="Maximum latitude in degrees for sample window")
    parser.add_argument("--sample-lon-min-deg", type=float, default=None,
                        help="Minimum longitude in degrees for sample window")
    parser.add_argument("--sample-lon-max-deg", type=float, default=None,
                        help="Maximum longitude in degrees for sample window")
    parser.add_argument("--sample-size", type=int, default=50,
                        help="Sample window size (default: 50)")
    parser.add_argument("--random-samples", type=int, default=0,
                        help="Number of random sample windows to include (default: 0)")
    parser.add_argument("--random-window", type=int, default=50,
                        help="Random sample window size (default: 50)")
    parser.add_argument("--random-seed", type=int, default=None,
                        help="Random seed for reproducible sampling")
    return parser.parse_args()


def sample_slice(
    var,
    lat_len,
    lon_len,
    time_len,
    lat_start=None,
    lon_start=None,
    size=50,
    lat_min=None,
    lat_max=None,
    lon_min=None,
    lon_max=None,
):
    # Default to mid-latitudes and first lon columns if not specified.
    t_end = min(time_len, 2)
    if lat_min is not None or lat_max is not None:
        lat_start = 0 if lat_min is None else max(0, lat_min)
        lat_end = lat_len if lat_max is None else min(lat_len, lat_max)
    else:
        if lat_start is None:
            lat_start = max(0, lat_len // 2 - size // 2)
        lat_end = min(lat_len, lat_start + size)

    if lon_min is not None or lon_max is not None:
        lon_start = 0 if lon_min is None else max(0, lon_min)
        lon_end = lon_len if lon_max is None else min(lon_len, lon_max)
    else:
        if lon_start is None:
            lon_start = 0
        lon_end = min(lon_len, lon_start + size)
    return var[:t_end, lat_start:lat_end, lon_start:lon_end], lat_start, lat_end, lon_start, lon_end


def main():
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output) if args.output else Path(DEFAULT_OUTPUT_DIR) / "inspect_smoothed.txt"

    with nc.Dataset(input_path, "r") as ds:
        dims = {k: len(v) for k, v in ds.dimensions.items()}
        variables = list(ds.variables.keys())
        lon = ds["longitude"][:]
        lat = ds["latitude"][:]
        time = ds["time"][:]

        report = []
        report.append(f"File: {input_path}")
        report.append(f"Format: {ds.file_format}")
        report.append(f"Dimensions: {dims}")
        report.append(f"Variables: {variables}")
        report.append(f"Longitude range: {float(lon.min())} to {float(lon.max())} (len={len(lon)})")
        report.append(f"Latitude range: {float(lat.min())} to {float(lat.max())} (len={len(lat)})")
        report.append(f"Time range: {float(time.min())} to {float(time.max())} (len={len(time)})")

        lulc = ds["LULC_states"]

        lat_min = args.sample_lat_min
        lat_max = args.sample_lat_max
        lon_min = args.sample_lon_min
        lon_max = args.sample_lon_max

        # Convert degree ranges to indices if provided
        if args.sample_lat_min_deg is not None or args.sample_lat_max_deg is not None:
            lat_min_deg = args.sample_lat_min_deg if args.sample_lat_min_deg is not None else lat.min()
            lat_max_deg = args.sample_lat_max_deg if args.sample_lat_max_deg is not None else lat.max()
            idx = np.where((lat >= lat_min_deg) & (lat <= lat_max_deg))[0]
            if len(idx) > 0:
                lat_min = int(idx.min())
                lat_max = int(idx.max()) + 1
        if args.sample_lon_min_deg is not None or args.sample_lon_max_deg is not None:
            lon_min_deg = args.sample_lon_min_deg if args.sample_lon_min_deg is not None else lon.min()
            lon_max_deg = args.sample_lon_max_deg if args.sample_lon_max_deg is not None else lon.max()
            idx = np.where((lon >= lon_min_deg) & (lon <= lon_max_deg))[0]
            if len(idx) > 0:
                lon_min = int(idx.min())
                lon_max = int(idx.max()) + 1

        sample, lat_start, lat_end, lon_start, lon_end = sample_slice(
            lulc, len(lat), len(lon), len(time),
            lat_start=args.sample_lat,
            lon_start=args.sample_lon,
            size=args.sample_size,
            lat_min=lat_min,
            lat_max=lat_max,
            lon_min=lon_min,
            lon_max=lon_max,
        )
        sample = np.asarray(sample)
        uniq = np.unique(sample)
        report.append(f"Sample window: lat[{lat_start}:{lat_end}] lon[{lon_start}:{lon_end}]")
        report.append(f"Sample dtype: {sample.dtype}")
        report.append(f"Sample min/max: {int(sample.min())} / {int(sample.max())}")
        report.append(f"Sample unique codes (first 20): {uniq[:20].tolist()}")
        report.append(f"Sample unique count: {len(uniq)}")

        if args.random_samples and args.random_samples > 0:
            if args.random_seed is not None:
                np.random.seed(args.random_seed)
            lat_len = len(lat)
            lon_len = len(lon)
            win = max(1, args.random_window)
            report.append(f"Random samples: count={args.random_samples}, window={win}")
            for i in range(args.random_samples):
                lat_start = np.random.randint(0, max(1, lat_len - win))
                lon_start = np.random.randint(0, max(1, lon_len - win))
                sample_r = np.asarray(lulc[:2, lat_start:lat_start+win, lon_start:lon_start+win])
                uniq_r = np.unique(sample_r)
                report.append(
                    f"Random {i+1}: lat[{lat_start}:{lat_start+win}] "
                    f"lon[{lon_start}:{lon_start+win}] "
                    f"unique_count={len(uniq_r)} min={int(uniq_r.min())} max={int(uniq_r.max())}"
                )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(report) + "\n")
    print(f"Wrote report to: {output_path}")


if __name__ == "__main__":
    main()
