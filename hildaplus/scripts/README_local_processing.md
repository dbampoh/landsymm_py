# Local Processing Guide for Updated HILDA+ Data

This guide explains how to process the updated HILDA+ data with new land cover categories on your local workstation.

## Quickstart (Minimal)

```bash
# Parallel full run (no FM, default output dir)
scripts/run_chain.sh --states /path/to/hildaplus_GLOB-2-1-crop_states.nc
```

## Overview

The processing involves two steps:
1. **Smoothing** the raw HILDA+ data to reduce flickering
2. **Upscaling** to generate LPJ-GUESS netfrac files at 0.5° resolution

## New Land Cover Categories in HILDA+ v2

The updated HILDA+ dataset includes these categories:
- 0: Ocean
- 11: Urban
- 22: Annual crops (NEW: separated from general cropland)
- 23: Tree crops (NEW)
- 24: Agroforestry (NEW)
- 33: Pasture/rangeland
- 40-45: Forest types
- 55: Unmanaged grass/shrubland
- 66: Sparse/no vegetation
- 77: Water
- 99: No data

Forest management codes (if in separate file):
- 0: No forest
- 1: Managed forest
- 2: Unmanaged forest

## Step 1: Smoothing Raw Data

### Script: `smoothing_local.py`

This script applies triple Gaussian filtering to reduce flickering in the time series.

### Usage:

```bash
# With forest management data in separate file
python smoothing_local.py --states hildap_states.nc --fm hildap_forest_mgmt.nc --output smoothed_output.nc

# With forest management integrated in states file (multiplied by 10)
python smoothing_local.py --states hildap_states.nc --output smoothed_output.nc

# Using default output name
python smoothing_local.py --states hildap_states.nc
```

#### Benchmark mode (estimate runtime)
Runs N chunks, reports speed, and estimates total runtime:
```bash
python smoothing_local.py --states hildap_states.nc --benchmark-chunks 3
```

Benchmark with a longitude range:
```bash
python smoothing_local.py --states hildap_states.nc \
  --benchmark-chunks 3 \
  --benchmark-lon-min -50 \
  --benchmark-lon-max -40
```

### Parallel option: `smoothing_local_parallel.py`
This version uses multiple CPU cores and writes safely to a single NetCDF.

```bash
python smoothing_local_parallel.py --states hildap_states.nc \
  --fm hildap_forest_mgmt.nc \
  --output /path/to/smoothed_output.nc \
  --workers 8 \
  --chunk-size 100
```

#### Benchmark mode (estimate runtime)
Runs N chunks, reports speed, and estimates total runtime:
```bash
python smoothing_local_parallel.py --states hildap_states.nc \
  --workers 8 \
  --chunk-size 100 \
  --benchmark-chunks 3
```

Benchmark with a longitude range:
```bash
python smoothing_local_parallel.py --states hildap_states.nc \
  --workers 8 \
  --chunk-size 100 \
  --benchmark-chunks 3 \
  --benchmark-lon-min -50 \
  --benchmark-lon-max -40
```

#### Create a matching mini-gridlist (for benchmark upscaling)
If you create a mini NetCDF from benchmark smoothing, filter your gridlist
to match its longitude/latitude extent:
```bash
python scripts/hilda_smoothing/make_minigridlist.py \
  --netcdf /path/to/mini_smoothed.nc \
  --gridlist /path/to/full_gridlist.txt \
  --output /path/to/mini_gridlist.txt
```

#### Inspect smoothed NetCDF
```bash
python scripts/hilda_smoothing/inspect_smoothed.py \
  --input /path/to/smoothed.nc
```

Inspect a specific lon/lat range (degrees):
```bash
python scripts/hilda_smoothing/inspect_smoothed.py \
  --input /path/to/smoothed.nc \
  --sample-lon-min-deg -10 --sample-lon-max-deg 10 \
  --sample-lat-min-deg 45 --sample-lat-max-deg 55
```

### Features:
- Memory-efficient chunked processing
- Progress reporting
- Handles both integrated and separate forest management data
- Optimized for single machine (no parallelization overhead)

### Performance Tips:
- Default chunk size is 1000 longitude pixels
- Adjust `CHUNK_SIZE` in the script if you have memory constraints
- Processing time depends on dataset size (expect ~1-2 hours for global data)

## Step 2: Generate Netfrac Files

### Script: `hildap_tables_netfrac_v3.py`

This updated script handles the new land cover categories and properly maps them to LPJ-GUESS categories.

### Land Cover Mapping:

| LPJ-GUESS Category | HILDA+ v2 Codes | Notes |
|-------------------|-----------------|--------|
| URBAN | 11 | Urban areas |
| CROPLAND | 22 + 23 + 24 | Annual crops + tree crops + agroforestry |
| PASTURE | 33 | Pasture/rangeland |
| FOREST | 40-45 + 400-450 | All forests (combined mode) |
| NATURAL | 55 + 66 | Grassland + sparse vegetation |
| BARREN | 0 + 77 + 99 | Ocean, water, no data |

**Important**:
- Default forest mode is **combined**, producing a single `FOREST` column.
- If you pass `--forest-mode split` and the smoothed file contains integrated FM codes (400–450), the output includes:
  `FOREST_MANAGED` + `FOREST_UNMANAGED` instead of `FOREST`.
- In split mode, `NATURAL` remains `55 + 66` only.

### Configuration:

You can either edit the script constants or pass paths via CLI flags.
CLI flags take precedence.
```bash
python hildap_tables_netfrac_v3.py \
  --datafile /path/to/smoothed_hildaplus.nc \
  --gridlist /path/to/gridlist.txt \
  --output /path/to/netfrac.txt \
  --forestfrac-output /path/to/forestfrac.txt \
  --forest-mode combined
```

### Usage:

```bash
# Process all grid points (using defaults in the script)
python hildap_tables_netfrac_v3.py

# Process subset (100 points starting at line 0)
python hildap_tables_netfrac_v3.py 0 100

# With output suffix
python hildap_tables_netfrac_v3.py 0 100 test
```

#### Benchmark mode (estimate runtime)
```bash
python hildap_tables_netfrac_v3.py --benchmark-lines 200
```

### Parallel option: `hildap_tables_netfrac_v3_parallel.py`
Splits the gridlist into chunks, runs workers in parallel, and merges outputs safely.

```bash
python hildap_tables_netfrac_v3_parallel.py \
  --workers 8 \
  --chunk-lines 500 \
  --forest-mode combined \
  --forestfrac-output /path/to/forestfrac.txt
```

#### Parallel benchmark mode
```bash
python hildap_tables_netfrac_v3_parallel.py \
  --workers 8 \
  --chunk-lines 500 \
  --benchmark-lines 200 \
  --forest-mode combined
```

#### Inspect netfrac output
```bash
python scripts/hildaplus-upscaling/inspect_netfrac_forestfrac.py \
  --input /path/to/netfrac.txt \
  --output /path/to/inspect_netfrac_forestfrac.txt \
  --gridlist /path/to/gridlist.txt
```

The inspection report includes netfrac sum checks, unique gridcell/year counts,
and forestfrac consistency checks (including `Max |FOREST*(sum forestfrac) - FOREST|`).

## Complete Workflow Example

```bash
# 1. Smooth the raw data
python smoothing_local.py raw_hildap_states.nc raw_hildap_fm.nc hildap_smoothed.nc

# 2. Create/obtain gridlist file
# Format: longitude latitude
# Example:
echo "-180.00 90.00" > my_gridlist.txt
echo "0.00 51.50" >> my_gridlist.txt
echo "-74.00 40.75" >> my_gridlist.txt

# 3. Update paths in hildap_tables_netfrac_v3.py

# 4. Generate netfrac file
python hildap_tables_netfrac_v3.py
```

## One-command Run Chain

Use the helper script to run smoothing + optional mini-gridlist + upscaling:

```bash
scripts/run_chain.sh --states /path/to/states.nc
```

Benchmark chain example:
```bash
scripts/run_chain.sh --states /path/to/states.nc --benchmark
```

Preview commands without running:
```bash
scripts/run_chain.sh --states /path/to/states.nc --benchmark --dry-run
```

Print resolved configuration only:
```bash
scripts/run_chain.sh --states /path/to/states.nc --print-config
```

## Full Run Examples (Technical How-To)

### A) Full run with individual scripts (single-process)
```bash
# 1) Smoothing (single-process, no FM)
python scripts/hilda_smoothing/smoothing_local.py \
  --states /path/to/hildaplus_GLOB-2-1-crop_states.nc \
  --output /path/to/output/hildaplus_smoothed.nc \
  --chunk-size 1000

# 2) Upscaling (single-process)
python scripts/hildaplus-upscaling/hildap_tables_netfrac_v3.py \
  --datafile /path/to/output/hildaplus_smoothed.nc \
  --gridlist /path/to/gridlist.txt \
  --output /path/to/output/hildaplus_netfrac_1901_2020.txt
```

### B) Full run with individual scripts (parallel)
```bash
# 1) Smoothing (parallel, no FM)
python scripts/hilda_smoothing/smoothing_local_parallel.py \
  --states /path/to/hildaplus_GLOB-2-1-crop_states.nc \
  --output /path/to/output/hildaplus_smoothed.nc \
  --workers 8 \
  --chunk-size 100

# 2) Upscaling (parallel)
python scripts/hildaplus-upscaling/hildap_tables_netfrac_v3_parallel.py \
  --datafile /path/to/output/hildaplus_smoothed.nc \
  --gridlist /path/to/gridlist.txt \
  --output /path/to/output/hildaplus_netfrac_1901_2020.txt \
  --workers 8 \
  --chunk-lines 500 \
  --forestfrac-output /path/to/output/hildaplus_forestfrac_1901_2020.txt \
  --forest-mode combined
### F) Full run with FM and split forests (parallel)
```bash
python scripts/hildaplus-upscaling/hildap_tables_netfrac_v3_parallel.py \
  --datafile /path/to/output/hildaplus_smoothed_with_fm.nc \
  --gridlist /path/to/gridlist.txt \
  --output /path/to/output/hildaplus_netfrac_1901_2020.txt \
  --workers 8 \
  --chunk-lines 500 \
  --forestfrac-output /path/to/output/hildaplus_forestfrac_1901_2020.txt \
  --forest-mode split
```
```

### C) Full run with FM file (single-process)
```bash
# 1) Smoothing (single-process, with FM)
python scripts/hilda_smoothing/smoothing_local.py \
  --states /path/to/hildaplus_GLOB-2-1-crop_states.nc \
  --fm /path/to/hildaplus_GLOB-2-0-f_forest-management.nc \
  --output /path/to/output/hildaplus_smoothed.nc \
  --chunk-size 1000

# 2) Upscaling (single-process, with FM)
python scripts/hildaplus-upscaling/hildap_tables_netfrac_v3.py \
  --datafile /path/to/output/hildaplus_smoothed.nc \
  --gridlist /path/to/gridlist.txt \
  --fmfile /path/to/hildaplus_GLOB-2-0-f_forest-management.nc \
  --output /path/to/output/hildaplus_netfrac_1901_2020.txt
```

### D) Full run via run-chain script (parallel)
```bash
scripts/run_chain.sh \
  --states /path/to/hildaplus_GLOB-2-1-crop_states.nc \
  --gridlist /path/to/gridlist.txt \
  --output-dir /path/to/output \
  --smooth-mode parallel \
  --upscale-mode parallel \
  --smooth-workers 8 \
  --smooth-chunk-size 100 \
  --upscale-workers 8 \
  --upscale-chunk-lines 500
```

### E) Full run via run-chain script (single-process)
```bash
scripts/run_chain.sh \
  --states /path/to/hildaplus_GLOB-2-1-crop_states.nc \
  --gridlist /path/to/gridlist.txt \
  --output-dir /path/to/output \
  --smooth-mode single \
  --upscale-mode single \
  --smooth-chunk-size 1000
```

## Output

The final output is `hildaplus_netfrac_1901_2020.txt` with columns:
- lon, lat, year
- URBAN, CROPLAND, PASTURE, FOREST, NATURAL, BARREN

Each row represents fractions (0.0-1.0) that should sum to ~1.0 (within rounding error).

## Troubleshooting

1. **Memory errors during smoothing**: Reduce `CHUNK_SIZE`
2. **Slow processing**: Normal for global datasets; consider processing regions separately
3. **Forest management issues**: Ensure FM data aligns temporally and spatially with states data
4. **Validation**: Check that fractions sum to 1.0 ± 0.0001 for each grid cell/year

## Validation Script

Here's a simple validation check:
```python
import numpy as np

# Load output
data = np.loadtxt('hildaplus_netfrac_1901_2020.txt', skiprows=1)

# Check sums (columns 3-8 are the fractions)
sums = data[:, 3:9].sum(axis=1)
bad_rows = np.where(np.abs(sums - 1.0) > 0.001)[0]

if len(bad_rows) > 0:
    print(f"Warning: {len(bad_rows)} rows don't sum to 1.0")
else:
    print("All rows sum to 1.0 ✓")
```