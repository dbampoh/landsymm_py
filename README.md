# LandSyMM Python Pipeline (`landsymm_py`)

Self-contained Python implementation of the **Land System Modular Model** pipeline,
covering all stages from raw HILDA+ processing through wetland/peatland integration.

## Pipeline Stages

| Stage | Package / Directory | Description |
|-------|---------------------|-------------|
| 1 | `hildaplus/` | HILDA+ smoothing, upscaling, and chain processing (standalone scripts) |
| 2 | `landsymm.remapping` | Baseline land-use construction from HILDA+ NetCDF → half-degree tables |
| 3 | `landsymm.harmonization` | PLUM scenario harmonization + conversion to LPJ-GUESS inputs |
| 4 | `landsymm.wetlands` | GLWD3 wetland/peatland integration into land-use data |

Shared utilities live in `landsymm.common`.

## Quick Start

### Installation

```bash
cd landsymm_py
pip install -e .
```

This installs the `landsymm` package in editable mode. All modules become importable
(e.g., `from landsymm.common.lpjg_io import read_table`).

### Data Directory

All data lives under `data/` at the project root:

```
landsymm_py/
├── data/
│   ├── geodata_py/          # Geospatial inputs
│   │   ├── HILDA+/data/     # HILDA+ input (raw NetCDF) + output (netfrac, forestfrac)
│   │   ├── MIRCA/           # MIRCA2000 crop fractions (originals/ + MIRCA.txt)
│   │   ├── AgGRID_*/        # N-fertilizer by crop (15 nc4 files)
│   │   ├── LUH2/            # LUH2 states (alternative to HILDA+)
│   │   ├── glwd3/           # GLWD3 wetland data + peatland_halfdeg.nc
│   │   ├── gridlists/       # Target gridlists
│   │   ├── soil/            # Soil input files (HWSD, soilmap)
│   │   └── country_boundaries/  # Country masks + codes (for figs)
│   ├── templates/                     # Excel templates for harm_by_numbers
│   ├── output_hildaplus_remap_10b/    # Remap baseline outputs
│   └── PLUMv2_LU_default_output/     # PLUM scenario data + harmonization outputs
```

Override the data directory by setting the `LANDSYMM_DATA_DIR` environment variable:

```bash
export LANDSYMM_DATA_DIR=/path/to/shared/data
```

### Running Each Stage

**Stage 1 — HILDA+ Smoothing & Upscaling:**

```bash
hildaplus/scripts/run_chain.sh --states /path/to/states.nc
```

**Stage 2 — Remapping:**

```bash
python -m landsymm.remapping.run_remap
```

**Stage 3 — PLUM Harmonization (full pipeline):**

```bash
python -m landsymm.harmonization.run_plumharm_pipeline
```

Or run individual steps:

```bash
python -m landsymm.harmonization.reformat_plum_gridded
python -m landsymm.harmonization.run_plumharm
python -m landsymm.harmonization.run_plumharm2lpjg
python -m landsymm.harmonization.run_plumharm_figs
```

**Stage 4 — Wetland/Peatland Integration:**

```bash
# Full pipeline (aggregation + insertion into HILDA+ remap LU)
python -m landsymm.wetlands.run_wetland_pipeline

# Or individual steps
python -m landsymm.wetlands.glwd3_to_halfdeg
python -m landsymm.wetlands.wetland_into_hilda
python -m landsymm.wetlands.wetland_into_forLPJG
```

All scripts accept `--help` for full option details.

## Project Structure

```
landsymm_py/
├── pyproject.toml                    # Package definition and dependencies
├── README.md                         # This file
├── landsymm_pipeline_technical_manual.md   # Comprehensive technical manual
├── forest_management_gross_transitions_inclusion_plan.md  # Future development plan
│
├── data/                             # All data (user-managed, not in version control)
│
├── landsymm/                         # Main Python package
│   ├── __init__.py
│   ├── config.py                     # Centralized path resolution
│   ├── common/                       # Shared utilities (I/O, interpolation, etc.)
│   ├── remapping/                    # Stage 2: baseline land-use construction
│   ├── harmonization/                # Stage 3: PLUM harmonization + LPJ-GUESS conversion
│   ├── wetlands/                     # Stage 4: GLWD3 wetland/peatland integration
│   └── tests/                        # Parity and validation tests
│
└── hildaplus/                        # Stage 1: HILDA+ processing (standalone scripts)
    └── scripts/
        ├── run_chain.sh              # Orchestrator for smoothing + upscaling
        ├── hilda_smoothing/          # Smoothing scripts
        └── hildaplus-upscaling/      # Upscaling scripts
```

## Dependencies

Core dependencies (installed automatically via `pip install -e .`):

- numpy, scipy, h5py, netCDF4, matplotlib, rasterio, xarray, pandas

Optional (for development):

```bash
pip install -e ".[dev]"   # adds pytest, ruff
```

## Configuration

Path resolution is centralized in `landsymm/config.py`. Key functions:

| Function | Returns |
|----------|---------|
| `get_project_root()` | `landsymm_py/` |
| `get_data_dir()` | `landsymm_py/data/` |
| `get_geodata_dir()` | `landsymm_py/data/geodata_py/` |
| `get_hildaplus_output_dir()` | `landsymm_py/data/geodata_py/HILDA+/data/output/` |
| `get_remap_output_dir()` | `landsymm_py/data/output_hildaplus_remap_10b/` |
| `get_plum_output_dir()` | `landsymm_py/data/PLUMv2_LU_default_output/` |

Usage in any module:

```python
from landsymm.config import get_data_dir, get_geodata_dir

data_path = get_data_dir() / "some_file.nc"
glwd3_dir = get_geodata_dir() / "glwd3"
```

## Documentation

- **Technical Manual:** `landsymm_pipeline_technical_manual.md` — comprehensive end-to-end documentation
- **Future Plans:** `forest_management_gross_transitions_inclusion_plan.md` — analysis of forest management and gross transitions integration
- **Per-stage READMEs:** in each subpackage directory (`remapping/`, `harmonization/`, `wetlands/`, `hildaplus/`, `common/`)
