#!/usr/bin/env python3
"""
Local version of HILDA+ smoothing script
Optimized for single machine processing
"""

import numpy as np
import netCDF4 as nc
from scipy import ndimage as nd
import argparse
import os
import sys
from pathlib import Path
from datetime import datetime

# Configuration
# Local-friendly defaults: sequential processing with configurable chunk sizes.
# The original script hardcodes cluster paths and writes per-slice .npz files.
SIGMA = 3.0  # Gaussian filter sigma value
CHUNK_SIZE = 1000  # Process this many longitude pixels at a time
OUTPUT_DIR = str(Path(__file__).resolve().parents[3] / 'data' / 'geodata_py' / 'HILDA+' / 'data' / 'output')  # Default output directory


def hildafomafusion(hildastates, hildafoma):
    """Merge HILDA+ LULC states with forest management data.
    Pixels where forest management is true are multiplied by 10.
    """
    return hildastates * (np.ones_like(hildafoma) + (hildafoma == 1) * 9)


def reduce_flickering_gaussian_1pixel(slice1, sigma=3.0):
    """Apply Gaussian filtering to reduce flickering in time series."""
    uniq = np.unique(slice1, return_inverse=True)
    if len(uniq[0]) > 1:  # Only process if more than 1 unique value
        # Create separate timeseries for all LULC classes
        u = [uniq[1] == i for i in range(len(uniq[0]))]
        # Apply Gaussian filtering
        gf = [nd.gaussian_filter(ui.astype(float), sigma, mode='constant', cval=0) for ui in u]
        # Reconstruct timeseries with maximum categories
        imax = np.argmax(gf, axis=0)
        newslice = slice1.copy()
        for i in range(len(uniq[0])):
            newslice[imax == i] = uniq[0][i]
        return newslice
    else:
        return slice1


def apply_gaussian_3times(slice1, sigma=3.0):
    """Apply Gaussian filtering three times (as per original)."""
    for _ in range(3):
        slice1 = reduce_flickering_gaussian_1pixel(slice1, sigma)
    return slice1


def smooth_hilda_data(hildastates_file, hildafoma_file=None, output_file=None,
                      chunk_size=CHUNK_SIZE, progress_interval=100, benchmark_chunks=0,
                      benchmark_lon_min=None, benchmark_lon_max=None):
    """
    Main smoothing function for local processing.
    
    Parameters:
    -----------
    hildastates_file : str
        Path to HILDA+ states NetCDF file
    hildafoma_file : str or None
        Path to forest management NetCDF file (optional)
    output_file : str
        Path for output smoothed NetCDF file
    chunk_size : int
        Number of longitude pixels to process at once
    progress_interval : int
        Report progress every N longitude chunks
    """
    
    print(f"Opening HILDA+ states file: {hildastates_file}")
    ds_states = nc.Dataset(hildastates_file, 'r')
    
    # Get dimensions
    time_dim = len(ds_states.dimensions['time'])
    lat_dim = len(ds_states.dimensions['latitude'])
    lon_dim = len(ds_states.dimensions['longitude'])
    
    print(f"Dataset dimensions: time={time_dim}, lat={lat_dim}, lon={lon_dim}")
    
    # Get coordinate variables
    times = ds_states.variables['time'][:]
    lats = ds_states.variables['latitude'][:]
    lons = ds_states.variables['longitude'][:]
    
    # Handle forest management data if provided (separate FM file).
    if hildafoma_file:
        print(f"Opening forest management file: {hildafoma_file}")
        ds_foma = nc.Dataset(hildafoma_file, 'r')
    else:
        print("No forest management file provided, proceeding without it")
        ds_foma = None
    
    # Resolve output path
    if output_file is None:
        output_file = os.path.join(OUTPUT_DIR, 'hildaplus_smoothed.nc')

    # Determine output longitude range for benchmark mode
    n_chunks = (lon_dim + chunk_size - 1) // chunk_size
    lon_start_idx = 0
    lon_end_idx = lon_dim
    if benchmark_chunks and benchmark_chunks > 0 and benchmark_lon_min is not None and benchmark_lon_max is not None:
        if benchmark_lon_min > benchmark_lon_max:
            benchmark_lon_min, benchmark_lon_max = benchmark_lon_max, benchmark_lon_min
        idx = np.where((lons >= benchmark_lon_min) & (lons <= benchmark_lon_max))[0]
        if len(idx) == 0:
            raise ValueError("Benchmark longitude range does not overlap dataset longitudes.")
        lon_start_idx = int(idx.min())
        lon_end_idx = int(idx.max()) + 1

    def build_chunks(start_idx, end_idx):
        chunk_start = start_idx // chunk_size
        chunk_end = (end_idx - 1) // chunk_size
        chunks = []
        for cidx in range(chunk_start, chunk_end + 1):
            lon_start = max(cidx * chunk_size, start_idx)
            lon_end = min((cidx + 1) * chunk_size, end_idx)
            chunks.append((cidx, lon_start, lon_end))
        return chunks

    chunks = build_chunks(lon_start_idx, lon_end_idx)
    if benchmark_chunks and benchmark_chunks > 0:
        chunks = chunks[:min(benchmark_chunks, len(chunks))]
        print(f"Benchmark mode: processing {len(chunks)} chunk(s)")

    out_lon_start = 0
    out_lon_end = lon_dim
    if benchmark_chunks and benchmark_chunks > 0:
        out_lon_start = chunks[0][1]
        out_lon_end = chunks[-1][2]
        print(f"Benchmark output will include longitude indices [{out_lon_start}:{out_lon_end}] only.")

    out_lon_dim = out_lon_end - out_lon_start

    # Create output NetCDF
    print(f"Creating output file: {output_file}")
    ds_out = nc.Dataset(output_file, 'w', format='NETCDF4')

    # Create dimensions
    ds_out.createDimension('time', time_dim)
    ds_out.createDimension('latitude', lat_dim)
    ds_out.createDimension('longitude', out_lon_dim)

    # Create coordinate variables
    time_var = ds_out.createVariable('time', 'i4', ('time',))
    lat_var = ds_out.createVariable('latitude', 'f4', ('latitude',))
    lon_var = ds_out.createVariable('longitude', 'f4', ('longitude',))

    # Copy coordinate data
    time_var[:] = times
    lat_var[:] = lats
    lon_var[:] = lons[out_lon_start:out_lon_end]
    
    # Copy attributes
    for attr in ds_states.variables['time'].ncattrs():
        if attr != '_FillValue':
            time_var.setncattr(attr, ds_states.variables['time'].getncattr(attr))
    for attr in ds_states.variables['latitude'].ncattrs():
        if attr != '_FillValue':
            lat_var.setncattr(attr, ds_states.variables['latitude'].getncattr(attr))
    for attr in ds_states.variables['longitude'].ncattrs():
        if attr != '_FillValue':
            lon_var.setncattr(attr, ds_states.variables['longitude'].getncattr(attr))
    
    # Create output variable
    lulc_var = ds_out.createVariable('LULC_states', 'u2', ('time', 'latitude', 'longitude'),
                                     zlib=True, complevel=4, chunksizes=(time_dim, 100, 100))
    lulc_var.setncattr('description', 'Smoothed HILDA+ land use/land cover states')
    lulc_var.setncattr('smoothing_method', 'Triple Gaussian filter with sigma=3.0')
    
    # Process in chunks to manage memory (single-machine friendly).
    print(f"\nProcessing {len(chunks)} chunk(s) of size {chunk_size}")
    start_time = datetime.now()

    for idx, (chunk_idx, lon_start, lon_end) in enumerate(chunks):
        if idx % progress_interval == 0:
            elapsed = (datetime.now() - start_time).total_seconds()
            percent = (idx / max(len(chunks), 1)) * 100
            print(f"Processing chunk {idx + 1}/{len(chunks)} ({percent:.1f}%) - Elapsed: {elapsed:.1f}s")
        
        # Read chunk data
        states_chunk = ds_states.variables['LULC_states'][:, :, lon_start:lon_end]
        
        if ds_foma:
            foma_chunk = ds_foma.variables['forest_management'][:, :, lon_start:lon_end]
            # Merge states with forest management
            data_chunk = hildafomafusion(states_chunk, foma_chunk)
        else:
            data_chunk = states_chunk
        
        # Apply smoothing to each pixel in the chunk (no parallelism here).
        smoothed_chunk = np.zeros_like(data_chunk)
        
        for i in range(data_chunk.shape[1]):  # latitude
            for j in range(data_chunk.shape[2]):  # longitude (within chunk)
                pixel_timeseries = data_chunk[:, i, j]
                smoothed_chunk[:, i, j] = apply_gaussian_3times(pixel_timeseries, SIGMA)
        
        # Write smoothed chunk to output (offset if benchmark range is used)
        out_start = lon_start - out_lon_start
        out_end = lon_end - out_lon_start
        lulc_var[:, :, out_start:out_end] = smoothed_chunk
    
    if benchmark_chunks and benchmark_chunks > 0:
        total_time = (datetime.now() - start_time).total_seconds()
        per_chunk = total_time / max(len(chunks), 1)
        est_total = per_chunk * n_chunks
        print(f"\nBenchmark complete: {len(chunks)} chunk(s) in {total_time:.1f}s "
              f"({per_chunk:.2f}s/chunk). Estimated full run: {est_total/3600:.2f} hours.")
        ds_states.close()
        if ds_foma:
            ds_foma.close()
        ds_out.close()
        return

    # Add global attributes
    ds_out.setncattr('title', 'Smoothed HILDA+ Land Use/Land Cover Dataset')
    ds_out.setncattr('source', f'Smoothed from: {os.path.basename(hildastates_file)}')
    ds_out.setncattr('history', f'Created on {datetime.now().isoformat()}')
    ds_out.setncattr('smoothing_sigma', SIGMA)
    
    # Close files
    ds_states.close()
    if ds_foma:
        ds_foma.close()
    ds_out.close()
    
    total_time = (datetime.now() - start_time).total_seconds()
    print(f"\nSmoothing completed in {total_time:.1f} seconds")
    print(f"Output saved to: {output_file}")


def parse_args():
    parser = argparse.ArgumentParser(description="Local HILDA+ smoothing (single process)")
    # Preferred flag-based args
    parser.add_argument("--states", default=None, help="Path to HILDA+ LULC states NetCDF")
    parser.add_argument("--fm", default=None, help="Path to forest management NetCDF (optional)")
    parser.add_argument("--output", default=None, help="Output NetCDF path (optional)")
    # Backward-compatible positional args
    parser.add_argument("states_pos", nargs="?", default=None,
                        help="(deprecated) Path to HILDA+ LULC states NetCDF")
    parser.add_argument("fm_pos", nargs="?", default=None,
                        help="(deprecated) Path to forest management NetCDF (optional)")
    parser.add_argument("output_pos", nargs="?", default=None,
                        help="(deprecated) Output NetCDF path (optional)")
    parser.add_argument("--chunk-size", type=int, default=CHUNK_SIZE,
                        help=f"Longitude chunk size (default {CHUNK_SIZE})")
    parser.add_argument("--benchmark-chunks", type=int, default=0,
                        help="Process N chunks, report timing, and exit (0 disables)")
    parser.add_argument("--benchmark-lon-min", type=float, default=None,
                        help="Benchmark longitude minimum (deg)")
    parser.add_argument("--benchmark-lon-max", type=float, default=None,
                        help="Benchmark longitude maximum (deg)")
    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_args()
    states = args.states or args.states_pos
    fm = args.fm if args.fm is not None else args.fm_pos
    output = args.output if args.output is not None else args.output_pos

    if not states:
        raise SystemExit("Missing required input: --states <path>")

    smooth_hilda_data(
        states,
        fm,
        output,
        chunk_size=args.chunk_size,
        benchmark_chunks=args.benchmark_chunks,
        benchmark_lon_min=args.benchmark_lon_min,
        benchmark_lon_max=args.benchmark_lon_max,
    )


if __name__ == '__main__':
    main()