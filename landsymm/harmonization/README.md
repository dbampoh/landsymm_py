# landsymm.harmonization

Python refactor of the MATLAB PLUM harmonization pipeline (`plum_harmonization/`).

This package harmonizes PLUM scenario land-use outputs against a historical baseline (HILDA+ remap) and converts them to LPJ-GUESS input format. It is a faithful port of the MATLAB codebase with validated parity.

---

## Pipeline Overview

The pipeline has three stages, matching the MATLAB workflow:

1. **PLUMharm** (`plumharm.py`) — Harmonizes raw PLUM scenario outputs against the HILDA+ baseline. Reads per-year PLUM inputs (LandCoverFract, CropFract, Fert, Irrig), applies area harmonization with ring redistribution, handles management (fertilization, irrigation), and produces per-year `.mat` output files.

2. **PLUMharm2LPJG** (`plumharm2lpjg.py`) — Converts the harmonized `.mat` outputs into LPJ-GUESS input text files (`landcover.txt`, `cropfractions.txt`, `nfert.txt`, `irrig.txt`).

3. **PLUMharmFigs** (`plumharm_figs.py`) — Generates diagnostic figures: time series, scatter plots, spatial maps, regional summaries, and optionally GeoTIFFs.

---

## Quick Start

Data locations and naming are resolved by `landsymm.config` and can be set with
environment variables (no source edits needed; see the Configuration section).
Scenarios default to every scenario directory found under the PLUM parent dir
(any subdirectory that contains a `{member}` dir, e.g. `s1/`); pass `--scenarios`
to restrict.

```bash
# Point the pipeline at your data / baseline / scenario set via env vars
export LANDSYMM_DATA_DIR=/path/to/data            # data root
export LANDSYMM_PLUM_DIRNAME=PLUMv2_LU_NFF_output  # PLUM scenario parent dir name
export LANDSYMM_REMAP_DIRNAME=output_hildaplus_remap_10b_3_py  # remap baseline dir name
# (optional) LANDSYMM_GEODATA_DIR, LANDSYMM_REMAP_VER, LANDSYMM_MEMBER
```

### Run the full pipeline (reformat + harmonization + conversion)

```bash
# All scenarios discovered under the parent dir, years 2021-2100
python -m landsymm.harmonization.run_plumharm_pipeline

# Single scenario, custom year range
python -m landsymm.harmonization.run_plumharm_pipeline --scenarios BAU --year1 2021 --yearN 2050

# Explicit parent dir (overrides LANDSYMM_PLUM_DIRNAME)
python -m landsymm.harmonization.run_plumharm_pipeline --parent-dir /path/to/PLUMv2_LU_NFF_output

# Skip reformatting (if already done)
python -m landsymm.harmonization.run_plumharm_pipeline --skip-reformat
```

### Run individual stages

```bash
# Harmonization only
python -m landsymm.harmonization.run_plumharm --scenarios BAU

# Conversion to LPJ-GUESS inputs only (requires harmonization outputs)
python -m landsymm.harmonization.run_plumharm2lpjg --scenarios BAU

# Diagnostic figures only (requires harmonization outputs)
python -m landsymm.harmonization.run_plumharm_figs --scenarios BAU
```

Note: remapping (Stage 1+2, `landsymm.remapping.run_remap`) and the optional
wetland/peatland step (Stage 4, `landsymm.wetlands.run_wetland_pipeline`) are
separate entrypoints and are NOT run by this harmonization pipeline.

### Scenarios

Scenarios are auto-discovered from the PLUM parent dir
(`config.discover_scenarios`), so any set works (the historical SSP set
`SSP1_RCP26` ... `SSP5_RCP85`, the NFF set `BAU`, `NfN_NfN`, ..., etc.).

---

## Input Data

### PLUM scenario inputs

Reformatted PLUM outputs, one directory per year:

```
data/PLUMv2_LU_default_output/{SSP}/s1/{YYYY}/
    LandCoverFract.txt    # Land cover fractions (NATURAL, CROPLAND, PASTURE, BARREN, URBAN)
    CropFract.txt         # Crop fractions within cropland
    Fert.txt              # Fertilization rates (kgN/ha)
    Irrig.txt             # Irrigation fractions
```

These are produced by `reformat_plum_gridded.py` (or its R predecessor `reformat_gridded_updated.R`) from raw PLUM outputs.

### Baseline reference data

Historical land-use baseline from the remapping pipeline. The directory name is
set by `LANDSYMM_REMAP_DIRNAME` and the version tag by `LANDSYMM_REMAP_VER`
(`remaps_v{ver}` subdir, `LU.remapv{ver}.txt` etc.); the helper
`config.get_remap_baseline_files()` returns these three paths:

```
<data>/<LANDSYMM_REMAP_DIRNAME>/remaps_v{ver}/
    LU.remapv{ver}.txt         # Land-use areas
    cropfracs.remapv{ver}.txt  # Crop fractions
    nfert.remapv{ver}.txt      # Fertilization rates
```

### Geodata

Static geographic datasets:

```
data/geodata_py/
    staticData_quarterdeg.nc    # Grid cell areas and land/water fractions
    maxcropfrac2.txt            # Terrain-limited maximum crop fraction
    protected_areas_with_points.txt  # Protected area fractions
    continents_from_countries/  # Continent shapefiles (for figures)
    country_boundaries/         # Country boundary data (for figures)
```

---

## Output Data

### PLUMharm outputs

Per-year `.mat` files with harmonized land-use and management data:

```
data/PLUMv2_LU_default_output/{SSP}/
    s1.HILDA+_remap_v10_old_62892_gL.harm.allow_unveg_py/{YYYY}/
        LandCoverFract.base2020.mat       # Half-degree land cover fractions
        LandCoverFract.base2020.2deg.mat  # 2-degree land cover fractions
        CropFract.base2020.mat            # Half-degree crop fractions
        CropFract.base2020.2deg.mat       # 2-degree crop fractions
        Fert.base2020.mat                 # Half-degree fertilization (kgN/ha)
        Fert.base2020.2deg.mat            # 2-degree fertilization
        Irrig.base2020.mat                # Half-degree irrigation fraction
        Irrig.base2020.2deg.mat           # 2-degree irrigation fraction
        post.base2020.mat                 # Carry state for restart
```

### PLUMharm2LPJG outputs

LPJ-GUESS input text files:

```
data/PLUMv2_LU_default_output/{SSP}/
    s1.HILDA+_remap_v10_old_62892_gL.harm.allow_unveg_py.forLPJG/
        landcover.txt        # Land cover fractions
        cropfractions.txt    # Crop fractions (irrigated + rainfed)
        nfert.txt            # Nitrogen fertilization (kgN/ha)
        irrig.txt            # Irrigation fractions
```

### PLUMharmFigs outputs

Diagnostic figures (PDF/PNG) and Excel summary tables:

```
data/PLUMv2_LU_default_output/{SSP}/
    s1.HILDA+_remap_v10_old_62892_gL.harm.allow_unveg_py_harms_figs/
        timeSeries_*.pdf       # Time series plots
        maps_*.png             # Spatial maps
        scatter_*.png          # Scatter plots
        harm_by_numbers.*.xlsx # Regional summary tables
```

Inline time-series figures (produced during harmonization):

```
    s1.HILDA+_remap_v10_old_62892_gL.harm.allow_unveg_py_figs/
        timeSeries_landUse.pdf
        timeSeries_crops.pdf
        timeSeries_nfert.pdf
        timeSeries_irrig.pdf
```

---

## Module Reference

### Core modules

| Module | MATLAB Equivalent | Description |
|---|---|---|
| `plumharm.py` | `PLUMharm.m` | Main harmonization loop |
| `plumharm_options.py` | `PLUMharm_options.m` | Configuration dataclass |
| `plumharm_import_ref.py` | `PLUMharm_importRefData.m` | Import baseline reference data |
| `plumharm_process_plum.py` | `PLUMharm_processPLUMin_areaCrops.m` | Process PLUM inputs |
| `plumharm_pp_read_plum.py` | `PLUMharm_pp_readPLUM.m` | Read/cache multi-year PLUM data |
| `plumharm_read_plum.py` | (wrapper) | Thin wrapper for pp_read_plum |
| `plumharm_area.py` | `PLUMharm_getUnmet_cropAreaRes.m`, `PLUMharm_fixTinyNegs.m` | Area unmet demand and tiny-neg fixes |
| `plumharm_dist.py` | `PLUMharm_distDeltas_areaCrops_recursive.m`, `PLUMharm_distMgmt.m`, `loop_thru_agri.m` | Delta distribution (2-deg to half-deg) |
| `plumharm_ring_redist.py` | `PLUMharm_ringRedist_areaCropsRes.m`, `PLUMharm_ringRedist_mgmt.m`, `PLUMharm_doRings_*.m` | Ring redistribution for area and management |
| `plumharm_mgmt.py` | `PLUMharm_getUnmet_mgmt.m`, `PLUMharm_interpolateMgmt.m` | Management harmonization (nfert, irrig) |
| `plumharm_checks.py` | `PLUMharm_checkCons_area.m`, `PLUMharm_checkCons_mgmt.m`, `PLUMharm_checkBadVals.m` | Conservation and sanity checks |
| `plumharm_debug.py` | `PLUMharm_debugOut.m`, `PLUMharm_debugOut_deltas.m` | Debug output helpers |

### Conversion and figures

| Module | MATLAB Equivalent | Description |
|---|---|---|
| `plumharm2lpjg.py` | `PLUMharm2LPJG.m` | Convert harmonized outputs to LPJ-GUESS inputs |
| `plumharm2lpjg_options.py` | `PLUMharm2LPJG_options.m` | Conversion configuration |
| `plumharm_figs.py` | `PLUMharmFigs.m`, `PLUMharmFigs_iterate_superReg.m`, `make_crops_timeseries_fig.m` | Diagnostic figure generation |
| `plumharm_figs_options.py` | `PLUMharmFigs_options.m` (not in repo) | Figures configuration |

### Run scripts

| Script | Description |
|---|---|
| `run_plumharm.py` | Run harmonization for all/selected scenarios |
| `run_plumharm2lpjg.py` | Run conversion for all/selected scenarios |
| `run_plumharm_figs.py` | Run diagnostic figures for all/selected scenarios |
| `run_plumharm_pipeline.py` | Run full pipeline (reformat + harmonization + conversion) |

### Shared dependencies (in `landsymm.common`)

| Module | Description |
|---|---|
| `lpjg_io.py` | LPJ-GUESS table I/O (read/write/maps2table) |
| `mapping_tools.py` | Grid mapping utilities |
| `interpolation.py` | Inpainting/interpolation (inpaint_nans) |
| `aggregation.py` | Grid aggregation (half-deg to 2-deg) |

---

## Dependencies

- Python 3.9+
- numpy
- scipy
- h5py
- netCDF4
- matplotlib
- pandas
- openpyxl (for Excel output in figures)
- rasterio (for GeoTIFF support in figures)
- geopandas (optional, for shapefile overlay in figures)

---

## Parity with MATLAB

This codebase has been validated against the MATLAB `plum_harmonization/` codebase:

- **LandCoverFract / CropFract**: Perfect parity (machine epsilon, ~1e-15)
- **Fert / Irrig**: 99.999%+ of values match. Remaining differences at 6 cells are bidirectional 1-ULP float64 arithmetic noise (4 of 6 are MATLAB artifacts, not Python)
- **PLUMharm2LPJG**: `landcover.txt` exact match; `cropfractions.txt` 8 last-digit-rounding diffs out of 106M values; `nfert.txt`/`irrig.txt` pass through upstream management noise only

See `plumharm_matlab_python_comprehensive_audit.md` in the project root for the full codebase comparison.

---

## Configuration

### Data locations and naming (environment variables)

Resolved by `landsymm.config` so no source edits are needed for a new data
location or scenario set:

| Env var | Default | Controls |
|---|---|---|
| `LANDSYMM_DATA_DIR` | `<project_root>/data` | data root |
| `LANDSYMM_GEODATA_DIR` | `<data>/geodata_py` | geodata dir (often elsewhere per machine) |
| `LANDSYMM_REMAP_DIRNAME` | `output_hildaplus_remap_10b` | remap-baseline dir name |
| `LANDSYMM_REMAP_VER` | `10_old_62892_gL` | remap version tag (baseline subdir + filenames) |
| `LANDSYMM_PLUM_DIRNAME` | `PLUMv2_LU_default_output` | PLUM scenario parent dir name |
| `LANDSYMM_MEMBER` | `s1` | ensemble member / median run |

Scenarios default to all discovered under the parent dir; override per run with
`--scenarios` / `--parent-dir`.

### Numeric / algorithmic options

These (matching MATLAB's `PLUMharm_options.m`) are set in the `run_*.py` scripts
via Python dataclasses:

| Parameter | Default | Description |
|---|---|---|
| `base_year` | 2020 | Baseline reference year |
| `year1` / `yearN` | 2021 / 2100 | Harmonization period |
| `allow_unveg` | True | Allow unvegetated cells |
| `norm2extra` | 0.177 | Fraction of crops allocated to ExtraCrop |
| `conserv_tol_pct` | 0.2 | Conservation tolerance (%) |
| `conserv_tol_area` | 1000 | Conservation tolerance (m²) |
| `fix_tiny_negs_tol_m2` | 1.0 | Tiny negative correction tolerance (m²) |
| `inpaint_method` | 4 | Inpainting method for missing management data |
| `out_prec` | 6 | Output decimal precision |
| `someofall` | True | Ensure minimum crop fraction everywhere |

To modify these numeric/algorithmic options, edit the relevant `run_*.py` script
directly. Data locations and the scenario set are controlled by the environment
variables above (and `--scenarios` / `--parent-dir`), not by source edits.
