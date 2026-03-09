#!/usr/bin/env python3
"""
Create a mini-gridlist that matches the lon/lat extent of a mini NetCDF output.
Intended for benchmark smoothing outputs with truncated longitude range.
"""

import argparse
from pathlib import Path

import netCDF4 as nc

DEFAULT_OUTPUT_DIR = str(Path(__file__).resolve().parents[3] / 'data' / 'geodata_py' / 'HILDA+' / 'data' / 'output')


def parse_args():
    parser = argparse.ArgumentParser(description="Filter gridlist to match NetCDF extent")
    parser.add_argument("--netcdf", required=True, help="Path to mini NetCDF (smoothed benchmark)")
    parser.add_argument("--gridlist", required=True, help="Path to full gridlist file")
    parser.add_argument(
        "--output",
        default=None,
        help=f"Path to write filtered mini-gridlist (default: {DEFAULT_OUTPUT_DIR}/mini_gridlist.txt)",
    )
    return parser.parse_args()


def read_bounds(nc_path):
    ds = nc.Dataset(nc_path, "r")
    lons = ds.variables["longitude"][:]
    lats = ds.variables["latitude"][:]
    ds.close()
    lon_min, lon_max = float(lons.min()), float(lons.max())
    lat_min, lat_max = float(lats.min()), float(lats.max())
    return lon_min, lon_max, lat_min, lat_max


def filter_gridlist(gridlist_path, output_path, lon_min, lon_max, lat_min, lat_max):
    total = 0
    kept = 0
    with open(gridlist_path, "r") as fin, open(output_path, "w") as fout:
        for line in fin:
            if not line.strip():
                continue
            total += 1
            parts = line.split()
            if len(parts) < 2:
                continue
            lon = float(parts[0])
            lat = float(parts[1])
            if lon_min <= lon <= lon_max and lat_min <= lat <= lat_max:
                fout.write(line)
                kept += 1
    return total, kept


def main():
    args = parse_args()
    netcdf_path = Path(args.netcdf)
    gridlist_path = Path(args.gridlist)
    output_path = Path(args.output) if args.output else Path(DEFAULT_OUTPUT_DIR) / "mini_gridlist.txt"

    lon_min, lon_max, lat_min, lat_max = read_bounds(netcdf_path)
    total, kept = filter_gridlist(gridlist_path, output_path, lon_min, lon_max, lat_min, lat_max)

    print(f"NetCDF bounds: lon [{lon_min}, {lon_max}], lat [{lat_min}, {lat_max}]")
    print(f"Gridlist lines: {total} total, {kept} kept")
    print(f"Mini-gridlist written to: {output_path}")


if __name__ == "__main__":
    main()
