# Wetlands-to-HILDA Pipeline (`landsymm.wetlands`)

Python pipeline for integrating GLWD3 wetland/peatland cover fractions into HILDA+ remap land-use data and PLUMharm2LPJG scenario output files. Ported from the R scripts in `wetlands_to_hilda/`.

---

## Overview

This package adds a **PEATLAND** land cover type to existing land-use datasets by carving a time-invariant peatland fraction from the **NATURAL** category, based on GLWD3 (Global Lakes and Wetlands Database, Level 3) high-resolution wetland data.

The package consists of three scripts that operate in sequence:

```
GLWD3 GeoTIFF (30 arc-second)
    │
    ▼  [Step 1: glwd3_to_halfdeg.py]
    │  Aggregate 12 GLWD3 classes to 0.5° fractions
    │  Apply Lehner (2004) scale factors
    │
    ▼
peatland_halfdeg.nc  (two wetland products: with/without forest types)
    │
    ├────────────────────────────────────────────┐
    │                                            │
    ▼  [Step 2a: wetland_into_hilda.py]          ▼  [Step 2b: wetland_into_forLPJG.py]
    │  Insert PEATLAND into HILDA+               │  Insert PEATLAND into scenario
    │  remap LU baseline file                    │  forLPJG landcover.txt files
    │                                            │
    ▼                                            ▼
LU.*_peatland.txt                        landcover_peatland.txt (per scenario)
```

### Approach H

All peatland insertion uses **Approach H** (from the original R scripts):

1. For each gridcell (Lon, Lat), find the **minimum NATURAL fraction** across all years.
2. Set `PEATLAND = min(min_NATURAL, GLWD3_wetland_frac)` — guaranteed to always be available from NATURAL.
3. `NATURAL_new = NATURAL_old - PEATLAND` at every year.
4. CROPLAND, PASTURE, BARREN are unchanged.
5. Fractions still sum to 1.0 at every gridcell and year.

This approach ensures PEATLAND never exceeds the available NATURAL, even at the year when NATURAL is at its minimum.

### Wetland Products

Two wetland products are available, selectable via `--wetland-product`:

| Product | Flag | Forest Classes | Use Case |
|---------|------|----------------|----------|
| **With forests** | `--wetland-product wforests` (default) | Includes Swamp Forest (class 5) and Coastal Wetland (class 6) | When vegetation dynamics model handles tree cover on wetlands internally |
| **Without forests** | `--wetland-product noforests` | Excludes classes 5 and 6 | When model explicitly represents swamp forests / mangroves in its PFTs |

---

## Prerequisites

### Data

- **GLWD3 GeoTIFF**: `data/geodata_py/glwd3/glwd_3.tif`
  - Source: [GLWD Level 3](https://www.worldwildlife.org/pages/global-lakes-and-wetlands-database) (Lehner & Döll 2004)
  - 30 arc-second resolution, 12 land-cover classes
- **HILDA+ remap LU file** (for Step 2a): `data/output_hildaplus_remap_10b/remaps_v10_old_62892_gL/LU.remapv10_old_62892_gL.txt`
- **Scenario forLPJG directories** (for Step 2b): `data/PLUMv2_LU_default_output/{SSP}/s1.*.forLPJG/landcover.txt`

### Python Packages

```
numpy
pandas
rasterio   (Step 1 only — reads GeoTIFF)
netCDF4    (Steps 1 & 2 — NetCDF I/O)
```

---

## Usage

### Step 1: Aggregate GLWD3 to Half-Degree

Reads the GLWD3 GeoTIFF and produces a NetCDF file with per-class fractions and total wetland fractions at 0.5° resolution.

```bash
python -m landsymm.wetlands.glwd3_to_halfdeg
```

**Options:**

| Flag | Default | Description |
|------|---------|-------------|
| `--glwd3-path` | `data/geodata_py/glwd3/glwd_3.tif` | Path to GLWD3 GeoTIFF |
| `--output-dir` | `data/geodata_py/glwd3/` | Output directory |
| `--quiet` | off | Suppress progress messages |

**Output:** `data/geodata_py/glwd3/peatland_halfdeg.nc`

Contains:
- `class_frac` — raw fraction of each 0.5° cell in each of 12 GLWD3 classes
- `scaled_frac_wforests` — scaled wetland fraction per class (with forest types)
- `scaled_frac_noforests` — scaled wetland fraction per class (without forest types)
- `wetland_frac` — total wetland fraction (with forest types)
- `wetland_frac_noforests` — total wetland fraction (without forest types)

### Step 2a: Insert Peatland into HILDA+ Remap LU

Produces a peatland-inclusive version of the HILDA+ remap LU baseline file.

```bash
python -m landsymm.wetlands.wetland_into_hilda
```

**Options:**

| Flag | Default | Description |
|------|---------|-------------|
| `--lu-path` | Auto-detected in `data/` | Path to HILDA+ remap LU file |
| `--wetland-nc` | `data/geodata_py/glwd3/peatland_halfdeg.nc` | Path to peatland NetCDF |
| `--wetland-product` | `wforests` | `wforests` or `noforests` |
| `--output-path` | Alongside input with `_peatland` tag | Output file path |
| `--quiet` | off | Suppress progress messages |

**Output:** `data/output_hildaplus_remap_10b/remaps_v10_old_62892_gL/LU.remapv10_old_62892_gL_peatland.txt`

Columns: `Lon Lat Year NATURAL CROPLAND PASTURE PEATLAND BARREN`

### Step 2b: Insert Peatland into Scenario forLPJG Landcover Files

Produces peatland-inclusive `landcover_peatland.txt` files alongside the original `landcover.txt` in each scenario's forLPJG directory.

```bash
# All scenarios
python -m landsymm.wetlands.wetland_into_forLPJG

# Specific scenarios
python -m landsymm.wetlands.wetland_into_forLPJG --scenarios SSP1_RCP26 SSP2_RCP45

# Custom parent directory
python -m landsymm.wetlands.wetland_into_forLPJG \
  --parent-dir data/PLUMv2_LU_default_output
```

**Options:**

| Flag | Default | Description |
|------|---------|-------------|
| `--parent-dir` | `data/PLUMv2_LU_default_output` | Parent directory with scenario folders |
| `--wetland-nc` | `data/geodata_py/glwd3/peatland_halfdeg.nc` | Path to peatland NetCDF |
| `--wetland-product` | `wforests` | `wforests` or `noforests` |
| `--scenarios` | all | Specific scenarios to process |
| `--quiet` | off | Suppress progress messages |

**Output:** `data/PLUMv2_LU_default_output/{SSP}/s1.*.forLPJG/landcover_peatland.txt`

Columns: `Lon Lat Year PASTURE CROPLAND NATURAL PEATLAND BARREN`

### Full Pipeline (Steps 1 + 2a)

```bash
python -m landsymm.wetlands.run_wetland_pipeline

# Skip Step 1 if peatland_halfdeg.nc already exists
python -m landsymm.wetlands.run_wetland_pipeline --skip-aggregation
```

---

## GLWD3 Class Mapping

| Class | Description | Scale Factor (with forests) | Scale Factor (without forests) |
|-------|-------------|:---------------------------:|:------------------------------:|
| 1 | Lake | 0.0 | 0.0 |
| 2 | Reservoir | 0.0 | 0.0 |
| 3 | River | 0.0 | 0.0 |
| 4 | Freshwater Marsh, Floodplain | 1.0 | 1.0 |
| 5 | Swamp Forest, Flooded Forest | 1.0 | **0.0** |
| 6 | Coastal Wetland | 1.0 | **0.0** |
| 7 | Pan, Brackish/Saline Wetland | 1.0 | 1.0 |
| 8 | Bog, Fen, Mire (Peatland) | 1.0 | 1.0 |
| 9 | Intermittent Wetland/Lake | 1.0 | 1.0 |
| 10 | 50-100% Wetland | 0.75 | 0.75 |
| 11 | 25-50% Wetland | 0.375 | 0.375 |
| 12 | Wetland Complex (0-25%) | 0.125 | 0.125 |

Scale factors are from Lehner & Döll (2004), Table 5.

---

## File Reference

| Script | R Original | Description |
|--------|------------|-------------|
| `glwd3_to_halfdeg.py` | `glwd3_to_halfdeg.R` | Aggregate GLWD3 to 0.5° wetland fractions |
| `wetland_into_hilda.py` | `glwd3_wetland_into_luh2.R` (Approach H) | Insert peatland into HILDA+ remap LU |
| `wetland_into_forLPJG.py` | — (new) | Insert peatland into scenario forLPJG landcover files |
| `run_wetland_pipeline.py` | — | Pipeline orchestration (Steps 1 + 2a) |

---

## Relationship to Original R Scripts

The original R scripts in `wetlands_to_hilda/` targeted LUH2 data and implemented multiple approaches (A-H) for wetland integration. This Python port:

- Targets **HILDA+ remap LU** data (and PLUMharm2LPJG scenario outputs) instead of LUH2
- Implements only **Approach H** (time-invariant peatland from minimum NATURAL)
- Adds the **forLPJG landcover.txt** processing step (not in original R scripts)
- Produces both "with forests" and "without forests" wetland products in a single run

---

## References

Lehner, B. & Döll, P. (2004). Development and validation of a global database of lakes, reservoirs and wetlands. *Journal of Hydrology*, 296, 1-22.
