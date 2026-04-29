# LandSyMM Python Pipeline (`landsymm_py`)

Self-contained Python implementation of the land-use pipeline for the
**Land System Modular Model (LandSyMM)**, covering all stages from raw
HILDA+ processing through baseline remapping, harmonization with PLUM
scenarios, and (optionally) wetland/peatland integration.

The pipeline constructs spatially explicit land-use forcing data for the
LPJ-GUESS dynamic global vegetation model. It bridges the gap between
historically observed land-use patterns (from the HILDA+ dataset at 1 km
resolution) and future land-use projections (from the PLUM land-use model)
by producing gridded inputs that transition seamlessly from the historical
period into scenario-driven futures. At each stage, multiple heterogeneous
data sources — HILDA+ land-use fractions, MIRCA2000 crop-specific harvested
areas, AgGRID nitrogen fertilization application rates, and (optionally)
GLWD3 wetland extents — are reconciled onto a common half-degree grid with
consistent land-use class definitions. The harmonization stage ensures that
the transition from observed baseline to projected scenario is spatially
conservative (total land area is preserved) and temporally smooth (no
discontinuities at the historical-to-future boundary), while respecting
constraints such as protected area fractions and minimum natural vegetation
rates. The final outputs are ready-to-use LPJ-GUESS input files — land cover
fractions, crop type allocations, nitrogen fertilization, and irrigation —
spanning both historical and future periods under multiple SSP-RCP scenarios.

## Why This Pipeline Exists — Scientific Context

`landsymm_py` is the data-preparation arm of the **Land System Modular
Model (LandSyMM)** framework, which couples a process-based dynamic global
vegetation model (LPJ-GUESS) with the Parsimonious Land Use Model
(PLUM/PLUMv2). The framework is described in two foundational publications:

- **Alexander et al. (2018)** — *Adaptation of global land use and management
  intensity to changes in climate and atmospheric carbon dioxide.* Global
  Change Biology 24:2791–2809. doi:10.1111/gcb.14110
- **Rabin et al. (2020)** — *Impacts of future agricultural change on
  ecosystem service indicators.* Earth System Dynamics 11:357–376.
  doi:10.5194/esd-11-357-2020

The LandSyMM workflow alternates LPJ-GUESS and PLUM runs to project how the
global agricultural system might adapt to climate, CO₂, and socioeconomic
change, and what the consequences are for ecosystem services such as carbon
storage, runoff, biodiversity, and nitrogen pollution. Within that workflow,
this pipeline is responsible for producing the **gridded land-use,
fertilizer, and irrigation forcing files that LPJ-GUESS reads** — both for
the historical period (1901–2014) and for SSP-RCP futures (2015–2100).

### Why each stage is needed

**Stage 1 (HILDA+ processing) and Stage 2 (Remapping):** LPJ-GUESS needs a
spatially complete, temporally consistent half-degree historical
land-use trajectory. HILDA+ provides observation-based land-cover fractions
at 1 km resolution; smoothing/upscaling and remapping convert these into
the LPJ-GUESS land-use class definitions (cropland, pasture, natural,
barren) at half-degree resolution. Crop-specific fractions from MIRCA2000
and nitrogen-fertilizer rates from AgGRID are reconciled onto the same
grid so that historical cropland is broken down by crop type and
fertilization intensity in a way that LPJ-GUESS can consume.

**Stage 3 (PLUM harmonization):** Quoting Rabin et al. (2020, Sect. 2.3):

> "The PLUM outputs must be processed first, because at the beginning of
> the future period they do not exactly match the land use and management
> forcings used at the end of the historical period. Feeding the raw PLUM
> outputs directly into LPJ-GUESS — causing large areas of sudden
> agricultural abandonment and expansion between 2010 and 2011 — would
> thus complicate interpretation of the results, especially of carbon
> cycling. We developed a harmonization routine, based on that published
> for Land Use Harmonization v1 dataset (LUH1) (Hurtt et al., 2011),
> that adjusts the PLUM outputs to ensure a smooth transition from the
> historical period to the future."

In other words, harmonization eliminates the spurious historical-to-future
discontinuity that would otherwise contaminate carbon-cycle attribution
in LandSyMM analyses. It is the necessary bridge between the
observation-derived baseline (Stage 2) and the model-projected future
(PLUM scenario outputs).

**Stage 4 (Wetland/peatland integration — optional):** This stage exists
specifically to support **coupled LPJ-GUESS ↔ IMOGEN simulations** in
which methane (CH₄) emissions from peatlands feed back into the IMOGEN
intermediate-complexity climate model alongside CO₂ and N₂O, and that
modified climate then drives subsequent LPJ-GUESS ecosystem responses.
GLWD3 wetland fractions are aggregated to half degree and inserted as a
distinct PEATLAND land-cover class in the otherwise non-wetland-aware
HILDA+ baseline (Stage 2 output) and in the harmonized PLUM scenario
output (Stage 3 output). **If you are not running coupled LPJ-GUESS ↔
IMOGEN simulations and do not need explicit peatland CH₄ accounting, you
can skip Stage 4 entirely** — the Stage 2 and Stage 3 outputs are
self-sufficient for vegetation/crop runs without peatland.

The integrated LPJ-GUESS LTS that consumes these outputs lives in a
companion repository:

- KIT GitLab: `https://gitlab.imk-ifu.kit.edu/bampoh-d/lpj-guess-integrated-landsymm`
- Helmholtz GitLab: `https://codebase.helmholtz.cloud/daniel.bampoh/lpj-guess-integrated-landsymm`
- GitHub: `https://github.com/dbampoh/LPJ-GUESS-integrated-LandSyMM`

When `landsymm_py` outputs are paired with that integrated LTS run with
its `iflandsymm_*` parameters enabled, the result reproduces the LandSyMM
fork's behavior end-to-end.

## Pipeline Stages

| Stage | Package / Directory | Required? | Purpose |
|-------|---------------------|-----------|---------|
| 1 | `hildaplus/` | Required | HILDA+ smoothing & upscaling — produces observation-based half-degree historical land-use trajectory (modified from Martin Winterbrink implementation) |
| 2 | `landsymm.remapping` | Required | Baseline LU construction from HILDA+ NetCDF + MIRCA crop fractions + AgGRID N-fert → half-degree LPJ-GUESS-class tables (modified from Sam Rabin LUH2 implementation) |
| 3 | `landsymm.harmonization` | Required for scenario runs | PLUM scenario harmonization with the historical baseline + conversion to LPJ-GUESS input format; eliminates the historical-to-future discontinuity that would corrupt carbon-cycle interpretation (modified from Sam Rabin LUH2 implementation) |
| 4 | `landsymm.wetlands` | **Optional** | GLWD3 peatland/wetland integration into LU files; needed only for coupled LPJ-GUESS ↔ IMOGEN simulations where peatland CH₄ emissions feed back into the climate model (modified from Peter Anthoni implementation) |

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

**Stage 4 — Wetland/Peatland Integration (optional):**

> **When to run this stage:** Only when you need explicit peatland
> representation in your LPJ-GUESS runs. This is typically required for
> **coupled LPJ-GUESS ↔ IMOGEN climate simulations** where peatland CH₄
> emissions feed back into the climate calculation alongside CO₂ and N₂O.
> If you are running standalone LPJ-GUESS with prescribed climate forcing
> and do not need explicit peatland CH₄ accounting, you can skip this
> stage — the Stage 2 / Stage 3 outputs are self-sufficient. The output
> files have a `_peatland` suffix to keep them separate from the
> non-peatland baseline.

```bash
# Full pipeline (aggregation + insertion into HILDA+ remap LU + scenario LU)
python -m landsymm.wetlands.run_wetland_pipeline

# Or individual steps
python -m landsymm.wetlands.glwd3_to_halfdeg       # GLWD3 1km/0.083° → 0.5°
python -m landsymm.wetlands.wetland_into_hilda     # Carve PEATLAND from NATURAL in Stage-2 baseline
python -m landsymm.wetlands.wetland_into_forLPJG   # Carve PEATLAND from NATURAL in each Stage-3 forLPJG/landcover.txt
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
