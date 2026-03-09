# LandSyMM Land-Use Pipeline — Technical Manual

**System**: HILDA+ Smoothing/Upscaling → Remapping → PLUM Harmonization → LPJ-GUESS Inputs → Wetland/Peatland Integration

---

## 1. System Overview

This system consists of three Python codebases that execute in succession to produce land-use inputs for LPJ-GUESS from raw HILDA+ satellite-derived land cover data and PLUM scenario projections.

```
┌──────────────────────────────────────────────────────────────────────┐
│  STAGE 1: hildaplus (smoothing & upscaling)                          │
│  Smoothing & Upscaling                                              │
│  Raw HILDA+ (0.01°) → Smoothed → Aggregated (0.5°)                 │
│  Output: hildaplus_netfrac_1901_2020.txt                            │
│          hildaplus_forestfrac_1901_2020.txt                          │
└────────────────────────┬─────────────────────────────────────────────┘
                         │
                         ▼
┌──────────────────────────────────────────────────────────────────────┐
│  STAGE 2: landsymm.remapping                                              │
│  Remapping & Baseline Construction                                  │
│  HILDA+ fractions + MIRCA crops + AgGRID nfert → Baseline LU       │
│  Output: LU.txt, cropfracs.txt, nfert.txt                           │
└────────────────────────┬─────────────────────────────────────────────┘
                         │
                         ▼
┌──────────────────────────────────────────────────────────────────────┐
│  PREPROCESSING: landsymm.harmonization/reformat_plum_gridded.py      │
│  Reformat raw PLUM outputs                                          │
│  LandCover.txt + LandUse.txt → LandCoverFract, CropFract, Fert,    │
│                                  Irrig (per year, per scenario)      │
└────────────────────────┬─────────────────────────────────────────────┘
                         │
                         ▼
┌──────────────────────────────────────────────────────────────────────┐
│  STAGE 3: landsymm.harmonization                                     │
│  PLUM Scenario Harmonization + LPJ-GUESS Conversion                 │
│  Baseline + PLUM scenarios → Harmonized outputs → forLPJG files     │
│  Output: landcover.txt, cropfractions.txt, nfert.txt, irrig.txt     │
└────────────────────────┬─────────────────────────────────────────────┘
                         │
                         ▼
┌──────────────────────────────────────────────────────────────────────┐
│  STAGE 4: landsymm.wetlands                                      │
│  Wetland/Peatland Integration                                       │
│  GLWD3 → 0.5° wetland fracs → PEATLAND into LU & forLPJG files    │
│  Output: LU.*_peatland.txt, landcover_peatland.txt (per scenario)   │
└──────────────────────────────────────────────────────────────────────┘
```

### Complete Execution Sequence

```bash
# Stage 1: Smoothing & Upscaling
cd hildaplus
bash scripts/run_chain.sh --states ... --gridlist ... --output-dir data/output

# Stage 2: Remapping
python -m landsymm.remapping.run_remap

# Stage 3: Reformat + Harmonization + LPJ-GUESS conversion (full pipeline)
python -m landsymm.harmonization.run_plumharm_pipeline \
  --parent-dir data/PLUMv2_LU_default_output

# Optional: Diagnostic figures
python -m landsymm.harmonization.run_plumharm_figs

# Stage 4: Wetland/Peatland Integration
# Step 1: Aggregate GLWD3 to half-degree
python -m landsymm.wetlands.glwd3_to_halfdeg

# Step 2a: Insert peatland into HILDA+ remap LU baseline
python -m landsymm.wetlands.wetland_into_hilda

# Step 2b: Insert peatland into all scenario forLPJG landcover files
python -m landsymm.wetlands.wetland_into_forLPJG \
  --parent-dir data/PLUMv2_LU_default_output
```

---

## 2. Stage 1: HILDA+ Smoothing & Upscaling

**Codebase**: `hildaplus/`

### Purpose

Processes raw HILDA+ land-use/land-cover (LULC) data at 0.01° resolution into aggregated 0.5° land-use fraction tables suitable for LPJ-GUESS.

### Pipeline

1. **Smoothing** — Applies triple Gaussian filtering (σ=3.0) to remove temporal "flickering" artifacts from the LULC classification time series. Each pixel's land cover type history is filtered independently.

2. **Upscaling** — Aggregates smoothed 0.01° data to 0.5° resolution by counting pixels per land cover category within each gridcell and computing fractions.

### Scripts

| Script | Purpose | Mode |
|---|---|---|
| `smoothing_local.py` | Smoothing (sequential) | Current |
| `smoothing_local_parallel.py` | Smoothing (multi-process) | Current |
| `hildap_tables_netfrac_v3.py` | Upscaling (sequential) | Current |
| `hildap_tables_netfrac_v3_parallel.py` | Upscaling (multi-process) | Current |
| `run_chain.sh` | Full pipeline orchestration | Current |
| `inspect_smoothed.py` | Validate smoothed NetCDF | Utility |
| `inspect_netfrac_forestfrac.py` | Validate upscaled tables | Utility |
| `make_minigridlist.py` | Create subset gridlists | Utility |
| `smoothingscript_run.py` | Smoothing (HPC cluster) | Legacy |
| `create_netcdf_from_slices.py` | Assemble HPC output slices | Legacy |
| `adapt_nc.py` | Add NetCDF metadata | Legacy |
| `hildap_tables.py` | Upscaling v1 (HPC) | Legacy |
| `launch_tables.py` | SLURM job generator | Legacy |

### Running

```bash
cd hildaplus
bash scripts/run_chain.sh \
  --states data/input_data/hildaplus_GLOB-2-1-crop_states.nc \
  --fm data/input_data/hildaplus_GLOB-2-0-f_forest-management.nc \
  --gridlist data/input_data/gridlist_in_62892_and_climate.txt \
  --output-dir data/output \
  --smooth-mode parallel \
  --upscale-mode parallel
```

### Inputs

| File | Description |
|---|---|
| `hildaplus_GLOB-2-1-crop_states.nc` | Raw HILDA+ LULC states (0.01°, 1900-2020) |
| `hildaplus_GLOB-2-0-f_forest-management.nc` | Forest management data (optional) |
| `gridlist_in_62892_and_climate.txt` | Target gridlist (62892 cells at 0.5°) |

### Outputs

Written to `data/geodata_py/HILDA+/data/output/` (configurable via `--output-dir`):

| File | Description |
|---|---|
| `hildaplus_netfrac_1901_2020.txt` | Net LU fractions per gridcell per year (URBAN, CROPLAND, PASTURE, FOREST, NATURAL, BARREN) |
| `hildaplus_forestfrac_1901_2020.txt` | Forest type sub-fractions (ForestNE, ForestND, ForestBE, ForestBD, ForestPNV) |
| `hildaplus_smoothed.nc` | Smoothed LULC NetCDF (intermediate) |

These outputs are consumed directly by Stage 2 (remapping) via `config.get_hildaplus_output_dir()`.

### Land Cover Mapping (v3, combined mode)

| Output Category | HILDA+ Codes |
|---|---|
| URBAN | 11 |
| CROPLAND | 22, 23 (tree crops), 24 (agroforestry) |
| PASTURE | 33 |
| FOREST | 40-45, 400-450 (managed) |
| NATURAL | 55, 66 |
| BARREN | 0, 77, 99 |

### Key Parameters

| Parameter | Default | Description |
|---|---|---|
| Gaussian σ | 3.0 | Smoothing kernel width |
| Repetitions | 3 | Number of smoothing passes |
| Target resolution | 0.5° | Output grid cell size |
| Precision | 7 decimals | Output text precision |
| Chunk size (serial) | 1000 lon pixels | Smoothing chunk |
| Workers (parallel) | 8 | Parallel processes |

### Dependencies

numpy, netCDF4, scipy

---

## 3. Stage 2: Remapping

**Codebase**: `landsymm.remapping/`

### Purpose

Constructs the historical baseline land-use dataset by combining HILDA+ land-use fractions with MIRCA crop allocation data and AgGRID fertilization data. Produces the reference inputs used by the harmonization stage.

### Pipeline

1. **Soil interpolation** — Fills missing values in soil input files using spring-analogy inpainting.
2. **Land-use import** — Reads HILDA+ (or LUH2) fractions, maps to 4 output types (NATURAL, CROPLAND, PASTURE, BARREN), adds ice/water fraction.
3. **Unvegetated handling** — Transfers tiny fraction from BARREN to NATURAL at fully barren cells.
4. **Sum normalization** — Forces LU fractions to sum exactly to 1.
5. **Crop fraction processing** — Reads MIRCA harvested areas, maps 26 crops to output crop list, computes ExtraCrop, normalizes.
6. **N-fertilizer** — Reads AgGRID data, maps to crops, converts units (kgN/ha → kgN/m²).
7. **Interpolation** — Inpaints missing cells for each crop using spring-analogy method.
8. **Output** — Writes LU, cropfracs, nfert tables and randomized gridlist.

### Running

```bash
python -m landsymm.remapping.run_remap
```

### Module Reference

| Module | MATLAB Equivalent | Description |
|---|---|---|
| `remap.py` | `remap.m` | Main pipeline driver |
| `remap_options.py` | `remap_options.m` | Configuration dataclass |
| `run_remap.py` | — | Entry point script |
| `lu_import.py` | `lu_import.m` | HILDA+/LUH2 import |
| `cropfrac.py` | `cropfrac.m` | MIRCA crop processing |
| `nfert.py` | `nfert.m` | AgGRID N-fertilizer |
| `soil.py` | `soil.m` | Soil interpolation |
| `diagnostics.py` | `diagnostics.m` | Missing-cell maps |

### Inputs

| File | Source | Description |
|---|---|---|
| `hildaplus_netfrac_1901_2020.txt` | Stage 1 (`geodata_py/HILDA+/data/output/`) | Land-use fractions |
| `MIRCA.txt` | Rebuilt from MIRCA2000 rasters (`geodata_py/MIRCA/.../originals/`) | Harvested area by crop (26 crops × RF/IR). A `.maps.mat` cache is used transparently for performance; if missing, `MIRCA.txt` is read directly. Can be regenerated via `python -m landsymm.remapping.rebuild_mirca_txt`. |
| `agmip_*_apprate_fill_NPK_0.5.nc4` | AgGRID | N-fertilizer by crop |
| `staticData_quarterdeg.nc` | External | Grid cell areas, land/water fractions |
| `gridlist_62892.*.txt` | External | Target output gridlist |
| Soil files (HWSD, soilmap) | External | Soil properties |

### Outputs

| File | Description |
|---|---|
| `LU.remapv{ver}.txt` | Land-use fractions (Lon, Lat, Year, NATURAL, CROPLAND, PASTURE, BARREN) |
| `cropfracs.remapv{ver}.txt` | Crop fractions per cell |
| `nfert.remapv{ver}.txt` | N-fertilizer rates per crop per cell (kgN/m²) |
| `gridlist.remapv{ver}.txt` | Output gridlist (randomized order) |
| Soil output files | Interpolated soil tables |

### Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `lu_source` | "HILDA+" | Land-use data source |
| `year_list_out` | 1901-2020 | Output year range |
| `inpaint_method` | 4 | Spring-analogy NaN filling |
| `plum_setaside_frac` | 0.088 | PLUM set-aside fraction |
| `fill_unveg` | 1e-6 | Fraction transferred to NATURAL at barren cells |
| `remap_ver` | "10_old_62892_gL" | Version tag for output filenames |

### Dependencies

numpy, scipy, netCDF4, h5py, pandas, matplotlib

---

## 4. PLUM Preprocessing (Reformatting)

**Script**: `landsymm.harmonization/reformat_plum_gridded.py`
**R original**: `data/PLUMv2_LU_default_output/reformat_gridded_updated.R`

### Purpose

Converts raw PLUM gridded outputs (LandCover.txt, LandUse.txt) into the harmonization-ready format expected by PLUMharm. This is a prerequisite step that runs before harmonization.

### What It Does

For each year directory under each scenario's `s1/` folder:

1. Reads `LandCover.txt` — extracts protected area fractions, aggregates land cover types (handling Photovoltaics → OtherNatural, Agrivoltaics → Cropland)
2. Reads `LandUse.txt` — normalizes per-crop management variables (FI, FQ, II, IQ, OI, Y) by area, converts crop areas to fractions of cropland
3. Renames original `LandUse.txt` to `LandUse_orig.txt` (preserves raw data)
4. Writes 5 output files:
   - `LandCoverFract.txt` — NATURAL, CROPLAND, PASTURE, BARREN, URBAN (fractions of total area)
   - `CropFract.txt` — per-crop fractions of cropland
   - `Fert.txt` — per-crop fertilization rates (FQ variable)
   - `Irrig.txt` — per-crop irrigation intensity (II variable)
   - `LandUse.txt` — reformatted combined table

### Crop Name Mapping

| PLUM Name | Output Name |
|---|---|
| wheat | CerealsC3 |
| oilcropsNFix | OilcropsNFix |
| oilcropsOther | OilcropsOther |
| starchyRoots | StarchyRoots |
| pulses | Pulses |
| maize | Maize |
| energycrops | Miscanthus |
| rice | Rice |
| fruitveg | FruitVeg |
| sugar | Sugar |
| setaside | Setaside |

### Running

```bash
# All scenarios
python -m landsymm.harmonization.reformat_plum_gridded \
  --parent-dir data/PLUMv2_LU_default_output

# Specific scenarios
python -m landsymm.harmonization.reformat_plum_gridded \
  --parent-dir data/PLUMv2_LU_default_output --scenarios SSP1_RCP26

# Single year directory
python -m landsymm.harmonization.reformat_plum_gridded \
  --single-year-dir data/PLUMv2_LU_default_output/SSP1_RCP26/s1/2021

# Custom parent directory
python -m landsymm.harmonization.reformat_plum_gridded \
  --parent-dir /path/to/PLUMv2_LU_default_output
```

### Input/Output

**Inputs** (per year directory):
```
{SSP}/{SSP}/s1/{YYYY}/
    LandCover.txt(.gz)    # Raw PLUM land cover (with Protection column)
    LandUse.txt(.gz)      # Raw PLUM land use (Crop × variable matrix)
```

**Outputs** (per year directory):
```
{SSP}/s1/{YYYY}/
    LandCoverFract.txt    # Land cover fractions
    CropFract.txt         # Crop fractions
    Fert.txt              # Fertilization rates
    Irrig.txt             # Irrigation intensity
    LandUse.txt           # Reformatted combined table
```

### Enhancements Over R Version

- CLI interface with `--scenarios`, `--single-year-dir`, `--quiet` flags
- Handles missing columns gracefully (Photovoltaics, Agrivoltaics may not exist in all PLUM versions)
- Gzip support for input files
- No R/dplyr/tidyr dependency

### Dependencies

numpy, pandas

---

## 5. Stage 3: PLUM Harmonization

**Codebase**: `landsymm.harmonization/`

### Purpose

Harmonizes PLUM scenario land-use projections (2021-2100) against the historical baseline from Stage 2. Converts harmonized outputs to LPJ-GUESS input format. Optionally generates diagnostic figures.

### Sub-stages

1. **PLUMharm** — Year-by-year harmonization of area and management at 2° resolution, distributed to 0.5°.
2. **PLUMharm2LPJG** — Stacks yearly `.mat` files into LPJ-GUESS text tables.
3. **PLUMharmFigs** — Diagnostic time series, maps, scatter plots, and regional summaries.

### Running

```bash
# Full pipeline (reformat + harmonization + conversion) — single scenario
python -m landsymm.harmonization.run_plumharm_pipeline \
  --parent-dir data/PLUMv2_LU_default_output --scenarios SSP1_RCP26

# Full pipeline — all scenarios
python -m landsymm.harmonization.run_plumharm_pipeline \
  --parent-dir data/PLUMv2_LU_default_output

# Skip reformatting if already done
python -m landsymm.harmonization.run_plumharm_pipeline \
  --parent-dir data/PLUMv2_LU_default_output --skip-reformat

# Individual stages (each also accepts --parent-dir)
python -m landsymm.harmonization.reformat_plum_gridded \
  --parent-dir data/PLUMv2_LU_default_output --scenarios SSP1_RCP26
python -m landsymm.harmonization.run_plumharm --scenarios SSP1_RCP26
python -m landsymm.harmonization.run_plumharm2lpjg --scenarios SSP1_RCP26
python -m landsymm.harmonization.run_plumharm_figs --scenarios SSP1_RCP26
```

### Module Reference

| Module | MATLAB Equivalent | Description |
|---|---|---|
| `plumharm.py` | `PLUMharm.m` | Main harmonization loop |
| `plumharm_options.py` | `PLUMharm_options.m` | Configuration |
| `plumharm_import_ref.py` | `PLUMharm_importRefData.m` | Baseline data import |
| `plumharm_process_plum.py` | `PLUMharm_processPLUMin_areaCrops.m` | PLUM input processing |
| `plumharm_area.py` | `PLUMharm_getUnmet_cropAreaRes.m` | Area unmet demand |
| `plumharm_dist.py` | `PLUMharm_distDeltas_*.m` | Delta distribution |
| `plumharm_ring_redist.py` | `PLUMharm_ringRedist_*.m` | Ring redistribution |
| `plumharm_mgmt.py` | `PLUMharm_getUnmet_mgmt.m` | Management harmonization |
| `plumharm_checks.py` | `PLUMharm_checkCons_*.m` | Conservation checks |
| `plumharm2lpjg.py` | `PLUMharm2LPJG.m` | LPJ-GUESS conversion |
| `plumharm_figs.py` | `PLUMharmFigs.m` | Diagnostic figures |

### Inputs

| File | Source | Description |
|---|---|---|
| `LU.remapv*.txt` | Stage 2 | Baseline land-use |
| `cropfracs.remapv*.txt` | Stage 2 | Baseline crop fractions |
| `nfert.remapv*.txt` | Stage 2 | Baseline fertilization |
| `{SSP}/s1/{YYYY}/*.txt` | PLUM (via R preprocessing) | Scenario LU/crop/fert/irrig |
| `staticData_quarterdeg.nc` | External | Grid areas |

### Outputs

| File | Description |
|---|---|
| `*.base2020.mat` / `*.2deg.mat` | Harmonized per-year land-use, crop, fert, irrig (half-deg and 2-deg) |
| `post.base2020.mat` | Restart state for year-over-year continuity |
| `landcover.txt` | LPJ-GUESS land cover input |
| `cropfractions.txt` | LPJ-GUESS crop fraction input (irrigated + rainfed) |
| `nfert.txt` | LPJ-GUESS N-fertilization input (kgN/ha) |
| `irrig.txt` | LPJ-GUESS irrigation input |
| `timeSeries_*.pdf` | Time series figures |

### Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `base_year` | 2020 | Baseline reference year |
| `year1` / `yearN` | 2021 / 2100 | Harmonization period |
| `allow_unveg` | True | Allow unvegetated cells |
| `norm2extra` | 0.177 | Fraction allocated to ExtraCrop |
| `conserv_tol_pct` | 0.2 | Conservation tolerance (%) |
| `fix_tiny_negs_tol_m2` | 1.0 | Tiny negative tolerance (m²) |
| `someofall` | True | Minimum crop fraction everywhere |
| `donation_order` | PASTURE, NATURAL, BARREN | Cropland donation order |

### Scenarios

| Code | Description |
|---|---|
| SSP1_RCP26 | Sustainability |
| SSP2_RCP45 | Middle of the road |
| SSP3_RCP70 | Regional rivalry |
| SSP4_RCP60 | Inequality |
| SSP5_RCP85 | Fossil-fueled development |

### Dependencies

numpy, scipy, h5py, netCDF4, matplotlib, pandas, openpyxl, rasterio (optional), geopandas (optional)

---

## 6. Stage 4: Wetland/Peatland Integration

**Codebase**: `landsymm.wetlands/`
**R originals**: `wetlands_to_hilda/glwd3_to_halfdeg.R`, `wetlands_to_hilda/glwd3_wetland_into_luh2.R`

### Purpose

Integrates GLWD3 (Global Lakes and Wetlands Database, Level 3) wetland/peatland cover fractions into the land-use datasets produced by Stages 2 and 3. Creates a new **PEATLAND** land cover type carved from the **NATURAL** category, producing peatland-inclusive versions of both the HILDA+ remap LU baseline file and each scenario's forLPJG landcover.txt.

### Pipeline

```
GLWD3 GeoTIFF (30 arc-second, 12 classes)
    │
    ▼  [Step 1: glwd3_to_halfdeg.py]
    │  Aggregate each class to 0.5° cell fractions
    │  Apply Lehner (2004) scale factors
    │  Produce "with forests" and "without forests" products
    │
    ▼
peatland_halfdeg.nc
    │
    ├───────────────────────────────────┐
    │                                   │
    ▼  [Step 2a: wetland_into_hilda.py] ▼  [Step 2b: wetland_into_forLPJG.py]
    │  HILDA+ remap LU baseline         │  Scenario forLPJG landcover files
    │  → LU.*_peatland.txt              │  → landcover_peatland.txt (×5 SSPs)
    │                                   │
    ▼                                   ▼
Peatland-inclusive LU baseline      Peatland-inclusive scenario inputs
(for future remapping runs)         (for LPJ-GUESS with peatland PFTs)
```

### Approach H — Time-Invariant Peatland from Minimum NATURAL

The peatland insertion algorithm (Approach H from the original R scripts):

1. For each gridcell `(Lon, Lat)`, compute `min_NATURAL` = minimum NATURAL fraction across **all** years in the dataset.
2. `PEATLAND = min(min_NATURAL, GLWD3_wetland_frac)` — guaranteed to never exceed available NATURAL at any year.
3. `NATURAL_new = NATURAL_old - PEATLAND` at every year.
4. CROPLAND, PASTURE, BARREN remain unchanged.
5. Row fractions still sum to 1.0 everywhere.

This guarantees non-negative NATURAL values for all years and gridcells.

### Wetland Products

Two wetland products are computed from the same GLWD3 source data:

| Product | CLI Flag | Classes Included | Use Case |
|---------|----------|------------------|----------|
| With forests | `--wetland-product wforests` (default) | All wetland classes (4-12) | Model handles tree cover on wetlands internally; forests are subsumed in NATURAL |
| Without forests | `--wetland-product noforests` | Classes 4, 7-12 (excludes Swamp Forest and Coastal Wetland) | Model has explicit swamp forest / mangrove PFTs that would double-count |

### GLWD3 Class Scale Factors

| Class | Description | Scale (with forests) | Scale (without forests) |
|-------|-------------|:--------------------:|:-----------------------:|
| 1 | Lake | 0.0 | 0.0 |
| 2 | Reservoir | 0.0 | 0.0 |
| 3 | River | 0.0 | 0.0 |
| 4 | Freshwater Marsh, Floodplain | 1.0 | 1.0 |
| 5 | Swamp Forest, Flooded Forest | 1.0 | 0.0 |
| 6 | Coastal Wetland | 1.0 | 0.0 |
| 7 | Pan, Brackish/Saline Wetland | 1.0 | 1.0 |
| 8 | Bog, Fen, Mire (Peatland) | 1.0 | 1.0 |
| 9 | Intermittent Wetland/Lake | 1.0 | 1.0 |
| 10 | 50-100% Wetland | 0.75 | 0.75 |
| 11 | 25-50% Wetland | 0.375 | 0.375 |
| 12 | Wetland Complex (0-25%) | 0.125 | 0.125 |

Scale factors from Lehner & Döll (2004), Table 5.

### Running

```bash
# Step 1: Aggregate GLWD3 to half-degree (produces peatland_halfdeg.nc)
python -m landsymm.wetlands.glwd3_to_halfdeg

# Step 2a: Insert peatland into HILDA+ remap LU baseline
python -m landsymm.wetlands.wetland_into_hilda

# Step 2a with "without forests" product
python -m landsymm.wetlands.wetland_into_hilda --wetland-product noforests

# Step 2b: Insert peatland into scenario forLPJG landcover files (all scenarios)
python -m landsymm.wetlands.wetland_into_forLPJG \
  --parent-dir data/PLUMv2_LU_default_output

# Step 2b: Specific scenarios only
python -m landsymm.wetlands.wetland_into_forLPJG \
  --parent-dir data/PLUMv2_LU_default_output --scenarios SSP1_RCP26 SSP3_RCP70

# Full pipeline (Steps 1 + 2a only)
python -m landsymm.wetlands.run_wetland_pipeline
python -m landsymm.wetlands.run_wetland_pipeline --skip-aggregation
```

### Module Reference

| Module | R Original | Description |
|--------|------------|-------------|
| `glwd3_to_halfdeg.py` | `glwd3_to_halfdeg.R` | Aggregate GLWD3 to 0.5° |
| `wetland_into_hilda.py` | `glwd3_wetland_into_luh2.R` (Approach H) | Peatland into HILDA+ remap LU |
| `wetland_into_forLPJG.py` | — (new) | Peatland into scenario forLPJG files |
| `run_wetland_pipeline.py` | — | Pipeline orchestration |

### Inputs

| File | Source | Description |
|------|--------|-------------|
| `glwd_3.tif` | GLWD3 (Lehner & Döll 2004) | 30 arc-second wetland classification |
| `LU.remapv*.txt` | Stage 2 | HILDA+ remap LU baseline |
| `landcover.txt` (per scenario) | Stage 3 | PLUMharm2LPJG scenario output |

### Outputs

| File | Location | Description |
|------|----------|-------------|
| `peatland_halfdeg.nc` | `data/geodata_py/glwd3/` | 0.5° wetland fractions (both products) |
| `LU.*_peatland.txt` | Alongside original LU file | Peatland-inclusive baseline LU |
| `landcover_peatland.txt` | Each scenario's forLPJG dir | Peatland-inclusive scenario landcover |

### Dependencies

numpy, pandas, rasterio, netCDF4

### Reference

Lehner, B. & Döll, P. (2004). Development and validation of a global database of lakes, reservoirs and wetlands. *Journal of Hydrology*, 296, 1-22.

---

## 7. Configuration & Path Resolution

**Module**: `landsymm/config.py`

All data paths are resolved relative to the project root. Override with the
`LANDSYMM_DATA_DIR` environment variable to use a custom data location.

### Path Functions

| Function | Default Path |
|----------|-------------|
| `get_project_root()` | `landsymm_py/` |
| `get_data_dir()` | `landsymm_py/data/` |
| `get_geodata_dir()` | `data/geodata_py/` |
| `get_hildaplus_output_dir()` | `data/geodata_py/HILDA+/data/output/` |
| `get_remap_output_dir()` | `data/output_hildaplus_remap_10b/` |
| `get_plum_output_dir()` | `data/PLUMv2_LU_default_output/` |

### Stage 1 → Stage 2 Data Flow

The HILDA+ smoothing/upscaling scripts (Stage 1) write their output to
`data/geodata_py/HILDA+/data/output/`. The remapping stage (Stage 2) reads
`hildaplus_netfrac_1901_2020.txt` from this same directory, creating a
direct pipeline connection:

```
hildaplus scripts → geodata_py/HILDA+/data/output/hildaplus_netfrac_1901_2020.txt
                                    ↓
landsymm.remapping.lu_import reads from geodata_py/HILDA+/data/output/
```

### Usage

```python
from landsymm.config import get_data_dir, get_geodata_dir, get_hildaplus_output_dir

data_path = get_data_dir() / "some_file.nc"
netfrac = get_hildaplus_output_dir() / "hildaplus_netfrac_1901_2020.txt"
```

---

## 8. Shared Utilities

**Codebase**: `landsymm/common/`

| Module | Lines | Description |
|---|---|---|
| `lpjg_io.py` | 1431 | LPJ-GUESS I/O: read/write text tables, MAT caching, maps↔table conversion |
| `interpolation.py` | 359 | NaN inpainting (6 methods, including spring-analogy PDE solver) |
| `mapping_tools.py` | 79 | Vector↔map conversions (Fortran-order for MATLAB parity) |
| `aggregation.py` | 36 | Resolution aggregation (block-sum) |
| `gridlist.py` | 40 | Gridlist loading |
| `io_utils.py` | 29 | Cache file cleanup |
| `geoarray.py` | 23 | GeoArray dataclass (stub) |
| `validation.py` | 16 | Sum-to-one assertion (stub) |
| `logging.py` | 12 | Logger setup (stub) |

---

## 9. Testing and Parity Validation

**Codebase**: `tests_py/`

| Script | Purpose | Usage |
|---|---|---|
| `plum_parity.py` | PLUMharm `.mat` file parity | `python -m tests_py.plum_parity --matlab-dir ... --python-dir ... --stats --max-cell --table output.md` |
| `plumharm2lpjg_parity.py` | PLUMharm2LPJG `.txt` file parity | `python -m tests_py.plumharm2lpjg_parity --matlab-dir ... --python-dir ... --stats --table output.md` |
| `run_parity.py` | Combined remap + PLUMharm parity | `python -m tests_py.run_parity --remap-matlab ... --remap-python ...` |
| `parity_utils.py` | Shared comparison utilities | (imported by other scripts) |

### Validated Parity Status

| Output | Status | Detail |
|---|---|---|
| LandCoverFract | **Perfect** | Machine epsilon (~1e-15) |
| CropFract | **Perfect** | ~1e-9 max (negligible) |
| Fert | **99.999%+** | 557/106M diffs from bidirectional 1-ULP noise |
| Irrig | **99.999%+** | 884/106M diffs from bidirectional 1-ULP noise |
| landcover.txt | **Perfect** | Exact match |
| cropfractions.txt | **Perfect** | 8 last-digit rounding diffs out of 106M |

---

## 10. End-to-End Data Flow

```
HILDA+ raw LULC (0.01°, 1900-2020)
    │
    ▼  [Stage 1: hildaplus]
    │  Smooth (3x Gaussian) → Upscale (0.01° → 0.5°)
    │
    ▼
hildaplus_netfrac_1901_2020.txt  (0.5°, 6 LU types, 120 years)
hildaplus_forestfrac_1901_2020.txt
    │
    ├──────────────────────────────┐
    ▼                              ▼
MIRCA harvested areas         AgGRID N-fertilizer
(MIRCA.txt, rebuilt from       (agmip_*_apprate_fill_NPK_0.5.nc4)
 52 .asc.gz rasters;            │
 .maps.mat cache used           │
 transparently)                 │
    │                              │
    ▼  [Stage 2: landsymm.remapping]     │
    │  Map crops, compute          │
    │  fractions, interpolate      │
    │  ◄──────────────────────────-┘
    │
    ▼
LU.remapv*.txt         (baseline LU fractions, 1901-2020)
cropfracs.remapv*.txt  (baseline crop fractions)
nfert.remapv*.txt      (baseline N-fertilization)
    │
    ├──────────────────────────────┐
    ▼                              ▼
Raw PLUM scenario outputs     staticData_quarterdeg.nc
(LandCover.txt, LandUse.txt)   (grid areas, land fractions)
    │
    ▼  [Preprocessing: reformat_plum_gridded.py]
    │  Reformat to LandCoverFract, CropFract, Fert, Irrig
    │
PLUM reformatted outputs      staticData_quarterdeg.nc
    │                              │
    ▼  [Stage 3: landsymm.harmonization]
    │  Harmonize areas (ring redistribution)
    │  Harmonize management (nfert, irrig)
    │  Convert to LPJ-GUESS format
    │
    ▼
landcover.txt          (LPJ-GUESS land cover input, 2021-2100)
cropfractions.txt      (LPJ-GUESS crop fractions)
nfert.txt              (LPJ-GUESS N-fertilization)
irrig.txt              (LPJ-GUESS irrigation)
    │
    ├──────────────────────────────┐
    ▼                              ▼
GLWD3 GeoTIFF (30 arc-sec)    LU.remapv*.txt (baseline)
    │                              │
    ▼  [Stage 4: landsymm.wetlands]
    │  Aggregate to 0.5° (Step 1)
    │  Insert PEATLAND (Approach H)
    │  into baseline (Step 2a) and
    │  forLPJG landcover (Step 2b)
    │
    ▼
landcover_peatland.txt (per scenario, with PEATLAND column)
LU.*_peatland.txt      (baseline with PEATLAND column)
    │
    ▼
LPJ-GUESS model runs (with or without peatland)
```

---

## 11. Directory Structure

```
landsymm_py/
├── pyproject.toml                      # Package definition & dependencies
├── README.md                           # Project overview
├── landsymm_pipeline_technical_manual.md  # This document
│
├── hildaplus/                          # Stage 1: Smoothing & upscaling (standalone)
│   ├── scripts/
│   │   ├── run_chain.sh                # Pipeline orchestration
│   │   ├── hilda_smoothing/            # Smoothing scripts
│   │   └── hildaplus-upscaling/        # Upscaling scripts
│   └── README.md                       # Documentation
│
├── landsymm/                           # Main Python package
│   ├── __init__.py
│   ├── config.py                       # Centralized path resolution
│   │
│   ├── common/                         # Shared utilities
│   │   ├── lpjg_io.py                  # LPJ-GUESS I/O
│   │   ├── interpolation.py            # NaN inpainting
│   │   ├── mapping_tools.py            # Grid conversions
│   │   └── aggregation.py              # Resolution aggregation
│   │
│   ├── remapping/                      # Stage 2: Remapping
│   │   ├── remap.py                    # Main driver
│   │   ├── run_remap.py                # Entry point
│   │   ├── rebuild_mirca_txt.py        # Reconstruct MIRCA.txt from raw rasters
│   │   └── README.md
│   │
│   ├── harmonization/                  # Stage 3: Harmonization
│   │   ├── plumharm.py                 # Main harmonization
│   │   ├── plumharm2lpjg.py            # LPJ-GUESS conversion
│   │   ├── reformat_plum_gridded.py    # PLUM preprocessing (R port)
│   │   ├── run_plumharm_pipeline.py    # Full pipeline
│   │   └── README.md
│   │
│   ├── wetlands/                       # Stage 4: Wetland integration
│   │   ├── glwd3_to_halfdeg.py         # GLWD3 aggregation
│   │   ├── wetland_into_hilda.py       # Peatland into baseline LU
│   │   ├── wetland_into_forLPJG.py     # Peatland into scenario landcover
│   │   ├── run_wetland_pipeline.py     # Pipeline orchestration
│   │   └── README.md
│   │
│   └── tests/                          # Parity and validation tests
│
└── data/                               # All input/output data
    ├── geodata_py/                     # Static geographic data
    │   ├── HILDA+/data/               # HILDA+ data
    │   │   ├── input_data/            # Raw HILDA+ NetCDF inputs
    │   │   └── output/                # Smoothed/upscaled outputs (netfrac, forestfrac)
    │   ├── MIRCA/harvested_area_grids_26crops_30mn/
    │   │   ├── originals/             # 52 raw .asc.gz rasters (source of truth)
    │   │   ├── MIRCA.txt              # Combined table (rebuilt from originals)
    │   │   └── MIRCA.txt.maps.mat     # Performance cache (auto-used by lpjg_io)
    │   ├── AgGRID_nutrient_input_v1.1/  # N-fertilizer nc4 files (15 crops)
    │   ├── LUH2/v2h/                 # LUH2 states (alternative to HILDA+)
    │   ├── gridlists/                 # Target gridlists
    │   ├── soil/                      # Soil input files (HWSD, soilmap)
    │   ├── country_boundaries/        # Country boundary masks + codes (for figs)
    │   ├── glwd3/                     # GLWD3 data + aggregated output
    │   ├── staticData_quarterdeg.nc   # Grid areas, land/water fractions
    │   ├── maxcropfrac2.txt           # Max crop fraction reference
    │   ├── protected_areas_with_points.txt  # Protected area mask
    │   └── wwf_terr_ecos_dissolveBiome_halfDeg_id.tif  # Biome mask (for figs)
    ├── templates/                      # Excel templates for harm_by_numbers
    ├── output_hildaplus_remap_10b/    # Stage 2 outputs (baseline)
    └── PLUMv2_LU_default_output/       # Stage 3 I/O (PLUM + harmonized)
```
