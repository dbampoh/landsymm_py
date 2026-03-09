#!/usr/bin/env python3
"""
Parallel local version of HILDA+ smoothing script.
Uses multiple CPU cores and a safe write strategy (temp chunks -> single NetCDF).
"""

import argparse
import os
from pathlib import Path
from datetime import datetime
from multiprocessing import Pool

import netCDF4 as nc
import numpy as np
from scipy import ndimage as nd

# Defaults
DEFAULT_SIGMA = 3.0
DEFAULT_CHUNK_SIZE = 100  # lon columns per chunk (matches original slice size)
DEFAULT_WORKERS = 8
DEFAULT_OUTPUT_DIR = str(Path(__file__).resolve().parents[3] / 'data' / 'geodata_py' / 'HILDA+' / 'data' / 'output')


def hildafomafusion(hildastates, hildafoma):
    """Merge HILDA+ LULC states with forest management data."""
    return hildastates * (np.ones_like(hildafoma) + (hildafoma == 1) * 9)


def reduce_flickering_gaussian_1pixel(slice1, sigma=3.0):
    """Apply Gaussian filtering to reduce flickering in a time series."""
    uniq = np.unique(slice1, return_inverse=True)
    if len(uniq[0]) > 1:
        u = [uniq[1] == i for i in range(len(uniq[0]))]
        gf = [nd.gaussian_filter(ui.astype(float), sigma, mode='constant', cval=0) for ui in u]
        imax = np.argmax(gf, axis=0)
        newslice = slice1.copy()
        for i in range(len(uniq[0])):
            newslice[imax == i] = uniq[0][i]
        return newslice
    return slice1


def apply_gaussian_3times(slice1, sigma=3.0):
    """Apply Gaussian filtering three times (as per original)."""
    for _ in range(3):
        slice1 = reduce_flickering_gaussian_1pixel(slice1, sigma)
    return slice1


def process_chunk(args):
    """Worker: smooth one longitude chunk and write to temp file."""
    (chunk_idx, lon_start, lon_end, states_path, fm_path, sigma, temp_dir) = args

    ds_states = nc.Dataset(states_path, 'r')
    states_chunk = np.array(ds_states['LULC_states'][:, :, lon_start:lon_end])

    if fm_path:
        ds_fm = nc.Dataset(fm_path, 'r')
        fm_chunk = np.array(ds_fm['forest_management'][:, :, lon_start:lon_end])
        data_chunk = hildafomafusion(states_chunk, fm_chunk)
        ds_fm.close()
    else:
        data_chunk = states_chunk

    smoothed_chunk = np.zeros_like(np.asarray(data_chunk))
    for i in range(data_chunk.shape[1]):  # latitude
        for j in range(data_chunk.shape[2]):  # longitude (within chunk)
            smoothed_chunk[:, i, j] = apply_gaussian_3times(data_chunk[:, i, j], sigma)

    ds_states.close()

    temp_path = os.path.join(temp_dir, f'chunk_{chunk_idx:04d}.npy')
    np.save(temp_path, smoothed_chunk)
    return lon_start, lon_end, temp_path


def parse_args():
    parser = argparse.ArgumentParser(description="Parallel HILDA+ smoothing (local)")
    # Preferred flag-based args
    parser.add_argument("--states", default=None, help="Path to HILDA+ LULC states NetCDF")
    parser.add_argument("--fm", default=None, help="Path to forest management NetCDF (optional)")
    parser.add_argument("--output", default=None, help="Output NetCDF path (optional)")
    # Backward-compatible positional args
    parser.add_argument("states_pos", nargs="?", default=None,
                        help="(deprecated) Path to HILDA+ LULC states NetCDF")
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE,
                        help=f"Longitude chunk size (default {DEFAULT_CHUNK_SIZE})")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                        help=f"Number of worker processes (default {DEFAULT_WORKERS})")
    parser.add_argument("--sigma", type=float, default=DEFAULT_SIGMA,
                        help=f"Gaussian sigma (default {DEFAULT_SIGMA})")
    parser.add_argument("--temp-dir", default=None,
                        help="Temp directory for chunk outputs (default: <output_dir>/tmp_chunks)")
    parser.add_argument("--benchmark-chunks", type=int, default=0,
                        help="Process N chunks, report timing, and exit (0 disables)")
    parser.add_argument("--benchmark-lon-min", type=float, default=None,
                        help="Benchmark longitude minimum (deg)")
    parser.add_argument("--benchmark-lon-max", type=float, default=None,
                        help="Benchmark longitude maximum (deg)")
    return parser.parse_args()


def main():
    args = parse_args()

    states_path = args.states or args.states_pos
    if not states_path:
        raise SystemExit("Missing required input: --states <path>")

    output_path = args.output or os.path.join(DEFAULT_OUTPUT_DIR, 'hildaplus_smoothed.nc')
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    temp_dir = args.temp_dir or os.path.join(os.path.dirname(output_path), 'tmp_chunks')
    os.makedirs(temp_dir, exist_ok=True)

    # Read metadata from states file
    ds_states = nc.Dataset(states_path, 'r')
    time_dim = len(ds_states.dimensions['time'])
    lat_dim = len(ds_states.dimensions['latitude'])
    lon_dim = len(ds_states.dimensions['longitude'])
    times = ds_states['time'][:]
    lats = ds_states['latitude'][:]
    lons = ds_states['longitude'][:]
    time_attrs = {a: ds_states['time'].getncattr(a) for a in ds_states['time'].ncattrs()}
    lat_attrs = {a: ds_states['latitude'].getncattr(a) for a in ds_states['latitude'].ncattrs()}
    lon_attrs = {a: ds_states['longitude'].getncattr(a) for a in ds_states['longitude'].ncattrs()}
    lulc_dtype = ds_states['LULC_states'].dtype
    ds_states.close()

    # Determine output longitude range for benchmark mode
    n_chunks = (lon_dim + args.chunk_size - 1) // args.chunk_size
    lon_start_idx = 0
    lon_end_idx = lon_dim
    if args.benchmark_chunks and args.benchmark_chunks > 0 and \
       args.benchmark_lon_min is not None and args.benchmark_lon_max is not None:
        lon_min = args.benchmark_lon_min
        lon_max = args.benchmark_lon_max
        if lon_min > lon_max:
            lon_min, lon_max = lon_max, lon_min
        idx = np.where((lons >= lon_min) & (lons <= lon_max))[0]
        if len(idx) == 0:
            raise ValueError("Benchmark longitude range does not overlap dataset longitudes.")
        lon_start_idx = int(idx.min())
        lon_end_idx = int(idx.max()) + 1

    def build_chunks(start_idx, end_idx):
        chunk_start = start_idx // args.chunk_size
        chunk_end = (end_idx - 1) // args.chunk_size
        chunks = []
        for cidx in range(chunk_start, chunk_end + 1):
            lon_start = max(cidx * args.chunk_size, start_idx)
            lon_end = min((cidx + 1) * args.chunk_size, end_idx)
            chunks.append((cidx, lon_start, lon_end, states_path, args.fm, args.sigma, temp_dir))
        return chunks

    chunks = build_chunks(lon_start_idx, lon_end_idx)
    if args.benchmark_chunks and args.benchmark_chunks > 0:
        chunks = chunks[:min(args.benchmark_chunks, len(chunks))]
        print(f"Benchmark mode: processing {len(chunks)} chunk(s) with {args.workers} workers")
    else:
        print(f"Processing {len(chunks)} chunk(s) with {args.workers} workers")

    out_lon_start = 0
    out_lon_end = lon_dim
    if args.benchmark_chunks and args.benchmark_chunks > 0:
        out_lon_start = chunks[0][1]
        out_lon_end = chunks[-1][2]
        print(f"Benchmark output will include longitude indices [{out_lon_start}:{out_lon_end}] only.")

    out_lon_dim = out_lon_end - out_lon_start

    # Prepare output NetCDF
    ds_out = nc.Dataset(output_path, 'w', format='NETCDF4')
    ds_out.createDimension('time', time_dim)
    ds_out.createDimension('latitude', lat_dim)
    ds_out.createDimension('longitude', out_lon_dim)
    time_var = ds_out.createVariable('time', 'i4', ('time',))
    lat_var = ds_out.createVariable('latitude', 'f4', ('latitude',))
    lon_var = ds_out.createVariable('longitude', 'f4', ('longitude',))
    time_var[:] = times
    lat_var[:] = lats
    lon_var[:] = lons[out_lon_start:out_lon_end]
    for k, v in time_attrs.items():
        if k != '_FillValue':
            time_var.setncattr(k, v)
    for k, v in lat_attrs.items():
        if k != '_FillValue':
            lat_var.setncattr(k, v)
    for k, v in lon_attrs.items():
        if k != '_FillValue':
            lon_var.setncattr(k, v)

    lulc_var = ds_out.createVariable('LULC_states', lulc_dtype,
                                     ('time', 'latitude', 'longitude'),
                                     zlib=True, complevel=4,
                                     chunksizes=(time_dim, 100, 100))
    lulc_var.setncattr('description', 'Smoothed HILDA+ land use/land cover states')
    lulc_var.setncattr('smoothing_method', f'Triple Gaussian filter sigma={args.sigma}')

    start_time = datetime.now()

    with Pool(processes=args.workers) as pool:
        for count, (lon_start, lon_end, temp_path) in enumerate(pool.imap_unordered(process_chunk, chunks), 1):
            data = np.load(temp_path)
            out_start = lon_start - out_lon_start
            out_end = lon_end - out_lon_start
            lulc_var[:, :, out_start:out_end] = data
            os.remove(temp_path)
            if count % 5 == 0 or count == len(chunks):
                elapsed = (datetime.now() - start_time).total_seconds()
                print(f"Written {count}/{len(chunks)} chunks - elapsed {elapsed:.1f}s")

    if args.benchmark_chunks and args.benchmark_chunks > 0:
        elapsed = (datetime.now() - start_time).total_seconds()
        per_chunk = elapsed / max(len(chunks), 1)
        est_total = per_chunk * n_chunks
        print(f"Benchmark complete: {len(chunks)} chunk(s) in {elapsed:.1f}s "
              f"({per_chunk:.2f}s/chunk). Estimated full run: {est_total/3600:.2f} hours.")
        ds_out.close()
        return

    ds_out.setncattr('title', 'Smoothed HILDA+ Land Use/Land Cover Dataset')
    ds_out.setncattr('source', f'Smoothed from: {os.path.basename(states_path)}')
    ds_out.setncattr('history', f'Created on {datetime.now().isoformat()}')
    ds_out.setncattr('smoothing_sigma', args.sigma)
    ds_out.close()

    total_time = (datetime.now() - start_time).total_seconds()
    print(f"Smoothing completed in {total_time:.1f} seconds")
    print(f"Output saved to: {output_path}")


if __name__ == '__main__':
    main()
