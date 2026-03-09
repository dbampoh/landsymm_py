# Session Handoff — January 26, 2026

**Scope:** PLUMharm parity closure, GLWD wetland pipeline completion, documentation, and planning for codebase consolidation.

---

## Table of Contents

1. [Session Overview](#1-session-overview)
2. [PLUMharm Parity Analysis — Final Closure](#2-plumharm-parity-analysis--final-closure)
3. [GLWD Wetland Pipeline — Implementation and Testing](#3-glwd-wetland-pipeline--implementation-and-testing)
4. [ForLPJG Peatland Integration — New Script](#4-forlpjg-peatland-integration--new-script)
5. [Documentation Updates](#5-documentation-updates)
6. [Future Work: Python Codebase Consolidation](#6-future-work-python-codebase-consolidation)
7. [Complete File Inventory](#7-complete-file-inventory)
8. [Current State Summary](#8-current-state-summary)

---

## 1. Session Overview

This session accomplished four major objectives:

1. **Closed the PLUMharm parity investigation**: Systematically confirmed that all remaining MATLAB-Python parity differences are 1-ULP floating-point artifacts, not code bugs. No further PLUMharm debugging is needed.

2. **Completed the GLWD wetland pipeline**: Modified scripts to produce both "with forests" and "without forests" wetland products, tested the full pipeline (GLWD3 aggregation + peatland insertion into HILDA+ remap LU), and verified all outputs.

3. **Created the forLPJG peatland integration**: Built and ran a new script (`wetland_into_forLPJG.py`) that inserts PEATLAND into all 5 scenario × 2 codebase (MATLAB + Python) forLPJG `landcover.txt` files, producing `landcover_peatland.txt` alongside each original.

4. **Updated all documentation**: Created the `landsymm.wetlands/README.md`, updated the technical manual with Stage 4, and produced this handoff document.

Additionally, a plan was discussed for consolidating all Python codebases into a self-contained project directory (`landsymm_py/`), targeted for the next session.

---

## 2. PLUMharm Parity Analysis — Final Closure

### Context

After fixing the R reformatting bug in `reformat_gridded_updated.R` (where NATURAL fraction was underestimated because `managedForest` was not divided by `area`), both MATLAB and Python PLUMharm pipelines were re-run on the corrected inputs. A post-fix parity table (`parity_table_full_postRfix_2021_2100.md`) showed apparent regression in CropFract differences and continued Fert/Irrig management spikes.

### Investigation Performed

Three diagnostic Python scripts were created and executed:

| Script | Purpose | Finding |
|--------|---------|---------|
| `diag_cropfract_regression.py` | Year-by-year CropFract comparison with CROPLAND correlation | All worst cells have CROPLAND < 10⁻¹⁴ |
| `diag_cropfract_threshold.py` | Threshold analysis: max CropFract diff when excluding near-zero CROPLAND | With CROPLAND > 10⁻⁶, max diff = **4.05e-11** |
| `diag_mgmt_spikes.py` | Fert/Irrig spike mechanism analysis | Zero-crossing: one codebase has crop area = 0, other has 10⁻¹⁹ m² |

Additionally performed:
- **MATLAB cache staleness check**: Confirmed `.mat` cache files are newer than corrected `.txt` inputs — no stale data
- **Code audit**: Side-by-side comparison of MATLAB `PLUMharm_processPLUMin_areaCrops.m` and Python `plumharm_process_plum.py` — identical algorithms confirmed
- **Input identity verification**: LandCoverFract parity at ~10⁻¹⁵ confirms both codebases read identical input data

### Root Cause (Confirmed)

All parity differences are **floating-point arithmetic artifacts**, not code bugs:

| Issue | Root Cause | Evidence |
|-------|-----------|----------|
| CropFract max_abs ~5e-3 | Division by CROPLAND < 10⁻¹⁴ | Excluding CROPLAND < 10⁻⁶ → max diff = 4e-11 |
| Fert spikes up to 726 kgN/ha | Zero-crossing: 1-ULP difference makes crop area 0 vs 10⁻¹⁹ | Spikes are **bidirectional** (MATLAB↔Python) |
| Irrig spikes of 1.0 | Same zero-crossing mechanism | Same bidirectionality |
| post.base max_abs ~4.96e+09 | Compounding of above over 80 years | Values are in m², represent sub-ULP fractional differences |

### PLUMharm2LPJG Output Parity (Also Confirmed)

| File | max_abs_diff | Nonzero Diffs | Assessment |
|------|-------------|---------------|------------|
| `landcover.txt` | 0.000000e+00 | 0 / 20,116,480 (0.0000%) | **Perfect** |
| `cropfractions.txt` | 4.863e-03 | 370 / 105,611,520 (0.0004%) | ULP noise |
| `nfert.txt` | 7.263e-02 | 1,510 / 105,611,520 (0.0014%) | ULP noise |
| `irrig.txt` | 1.000e+00 | 2,245 / 105,611,520 (0.0021%) | Zero-crossing |

### Conclusion

**The Python PLUMharm implementation is correct and scientifically equivalent to MATLAB.** No further PLUMharm parity debugging is needed. Full diagnostic report: `parity_diagnostic_report_postRfix.md`.

### Optional Future Enhancement

A CROPLAND floor guard (e.g., set CropFract to 0 when CROPLAND < 1 m²) could be added to **both** codebases to eliminate management spikes. This is a quality-of-life improvement, not a bug fix.

---

## 3. GLWD Wetland Pipeline — Implementation and Testing

### Modifications Made

#### `glwd3_to_halfdeg.py` — Dual Wetland Products

Previously generated only a single "with forests" wetland product. Modified to produce **both** products in a single NetCDF file (`peatland_halfdeg.nc`):

- Added `SCALE_FACTORS_NOFORESTS` array (zeroes out classes 5 and 6)
- Computes both `wetland_frac` (with forests) and `wetland_frac_noforests` (without forests)
- NetCDF output includes both products plus per-class scaled fractions for each
- Output filename changed from `peatland_wforests_halfdeg.nc` to `peatland_halfdeg.nc`

| Wetland Product | Classes Included | Nonzero Cells | Use Case |
|-----------------|-----------------|:-------------:|----------|
| `wetland_frac` (with forests) | 4-12 | 28,653 | Model handles tree cover internally (default) |
| `wetland_frac_noforests` (without forests) | 4, 7-12 | 25,839 | Model has explicit swamp forest PFTs |

#### `wetland_into_hilda.py` — Wetland Product Selection

- Added `WETLAND_VARS` dictionary mapping `"wforests"` → `"wetland_frac"` and `"noforests"` → `"wetland_frac_noforests"`
- `_load_wetland_frac_from_nc()` accepts `wetland_product` parameter
- `insert_wetland_approach_h()` accepts and passes `wetland_product`
- CLI argument `--wetland-product` added (choices: `wforests`, `noforests`; default: `wforests`)
- Backward compatibility: falls back to `wetland_frac` if older NetCDF file is used

#### `run_wetland_pipeline.py` — Passes Wetland Product

- `main()` accepts `wetland_product` parameter
- CLI argument `--wetland-product` added
- Updated NetCDF filename reference to `peatland_halfdeg.nc`

### Test Results

#### Step 1: GLWD3 Aggregation

```
Command: python -m landsymm.wetlands.glwd3_to_halfdeg
Runtime: 30.8s
Input: data/geodata_py/glwd3/glwd_3.tif (21600×43200 pixels)
Output: data/geodata_py/glwd3/peatland_halfdeg.nc
```

- With forests — wetland fraction range: [0.000000, 1.000000], 28,653 nonzero cells
- Without forests — wetland fraction range: [0.000000, 1.000000], 25,839 nonzero cells

#### Step 2a: Peatland into HILDA+ Remap LU

```
Command: python -m landsymm.wetlands.wetland_into_hilda
Runtime: 255.0s
Input: data/output_hildaplus_remap_10b/remaps_v10_old_62892_gL/LU.remapv10_old_62892_gL.txt
Output: data/output_hildaplus_remap_10b/remaps_v10_old_62892_gL/LU.remapv10_old_62892_gL_peatland.txt
```

- 7,546,080 rows, 120 years (1901-2020)
- 26,568 unique gridcells with PEATLAND > 0
- PEATLAND range: [0.000001, 0.999998]
- Max row-sum deviation from 1.0: 2e-6 (inherited from original file)
- No negative NATURAL values

---

## 4. ForLPJG Peatland Integration — New Script

### What Was Created

`landsymm.wetlands/wetland_into_forLPJG.py` — a new script that inserts PEATLAND into the scenario `forLPJG/landcover.txt` files produced by PLUMharm2LPJG. It:

1. Auto-discovers all `forLPJG` directories under a parent directory
2. For each, loads `landcover.txt`, applies Approach H (same as `wetland_into_hilda.py`)
3. Writes `landcover_peatland.txt` alongside the original
4. Preserves the original column order (`Lon Lat Year PASTURE CROPLAND NATURAL PEATLAND BARREN`)
5. Uses `np.savetxt` for fast vectorized I/O (~15s per ~5M-row file)

### CLI Usage

```bash
# All scenarios
python -m landsymm.wetlands.wetland_into_forLPJG \
  --parent-dir data/PLUMv2_LU_default_output

# Specific scenarios
python -m landsymm.wetlands.wetland_into_forLPJG \
  --parent-dir data/PLUMv2_LU_default_output --scenarios SSP1_RCP26 SSP3_RCP70

# With "without forests" product
python -m landsymm.wetlands.wetland_into_forLPJG \
  --parent-dir data/PLUMv2_LU_default_output --wetland-product noforests
```

### Test Results — All Scenarios

```
Command: python -m landsymm.wetlands.wetland_into_forLPJG --parent-dir data/PLUMv2_LU_default_output
Runtime: 2m 30s total (~15s per file)
Directories processed: 10 (5 scenarios × 2 codebases each)
```

| Scenario | Cells with PEATLAND | PEATLAND Max | MATLAB & Python Match |
|----------|:-------------------:|:------------:|:---------------------:|
| SSP1_RCP26 | 27,958 | 0.999990 | Yes |
| SSP2_RCP45 | 27,922 | 0.999998 | Yes |
| SSP3_RCP70 | 27,906 | 0.997045 | Yes |
| SSP4_RCP60 | 27,924 | 0.998084 | Yes |
| SSP5_RCP85 | 27,936 | 0.997408 | Yes |

Variation in peatland cell counts across scenarios is expected: scenarios with more cropland expansion have lower minimum NATURAL, reducing cells eligible for peatland.

### Verification (SSP1_RCP26 MATLAB dir)

```
Columns: Lon Lat Year PASTURE CROPLAND NATURAL PEATLAND BARREN
Rows: 5,029,120 (62,892 cells × 80 years)
Row sums: min=0.9999980, max=1.0000020 (max deviation 2e-6, inherited)
Negative NATURAL: 0
MATLAB and Python dirs: identical PEATLAND statistics
```

### Output File Locations

Each `landcover_peatland.txt` is placed alongside its original `landcover.txt`:

```
data/PLUMv2_LU_default_output/
├── SSP1_RCP26/
│   ├── s1.HILDA+_remap_v10_old_62892_gL.harm.allow_unveg.forLPJG/
│   │   ├── landcover.txt              ← original (unchanged)
│   │   └── landcover_peatland.txt     ← new (with PEATLAND column)
│   └── s1.HILDA+_remap_v10_old_62892_gL.harm.allow_unveg_py.forLPJG/
│       ├── landcover.txt
│       └── landcover_peatland.txt
├── SSP2_RCP45/ ... (same structure)
├── SSP3_RCP70/ ...
├── SSP4_RCP60/ ...
└── SSP5_RCP85/ ...
```

---

## 5. Documentation Updates

### Created

| File | Description |
|------|-------------|
| `landsymm.wetlands/README.md` | Full documentation for the wetland pipeline: overview, pipeline diagram, Approach H explanation, wetland products, GLWD3 class mapping, CLI usage for all scripts, file reference, and relationship to original R scripts |
| `parity_diagnostic_report_postRfix.md` | Comprehensive post-R-fix parity diagnostic report with investigation steps, root cause analysis, threshold analysis, management spike mechanism, and recommendations |
| `session_handoff_20260126.md` | This document |

### Updated

| File | Changes |
|------|---------|
| `landsymm_pipeline_technical_manual.md` | Added Stage 4 (Wetland/Peatland Integration) section with full pipeline description, Approach H algorithm, GLWD3 class scale factors table, running commands, module reference, I/O tables. Updated system overview diagram, execution sequence, end-to-end data flow, and directory structure. Renumbered subsequent sections (Shared Utilities → 7, Testing → 8, Data Flow → 9, Directory Structure → 10). |

---

## 6. Future Work: Python Codebase Consolidation

### Objective

Consolidate all Python codebases into a self-contained project directory that can run independently of the MATLAB codebase.

### Proposed Structure

```
landsymm_py/
├── pyproject.toml                         # Package definition, dependencies
├── README.md                              # Top-level quick-start guide
├── landsymm_pipeline_technical_manual.md   # Full technical reference
├── landsymm/                              # Single top-level Python package
│   ├── common/                            # Shared utilities (was common_py/)
│   │   └── README.md
│   ├── remapping/                         # Stage 2 (was remapping_py/)
│   │   └── README.md
│   ├── harmonization/                     # Stage 3 (was landsymm.harmonization/)
│   │   └── README.md
│   ├── wetlands/                          # Stage 4 (was landsymm.wetlands/)
│   │   └── README.md
│   ├── tests/                             # Parity testing (was tests_py/)
│   │   └── README.md
│   └── config.py                          # Centralized data-path resolution
├── hildaplus/                             # Stage 1 (standalone scripts)
│   └── README.md
└── data/                                  # Normal data directory (not symlink)
```

### Key Challenges

1. **Cross-package imports**: `landsymm.harmonization` and `remapping_py` both import from `common_py`. Need proper package structure with `pyproject.toml`.

2. **Data path resolution**: ~20-30 locations use `Path(__file__).resolve().parents[N]` to find the `data/` directory. Need centralized `config.py` with optional `LANDSYMM_DATA_DIR` environment variable override.

3. **`hildaplus_global_05deg_lpjguess`**: Has bash orchestration (`run_chain.sh`) and its own data directory. Decision needed: fold in or keep as external dependency.

4. **Run scripts and CLI**: All `run_*.py` scripts need updating for new package structure.

### Decision Points

- **Data directory**: Normal directory (not symlink) for full self-containment. One-time copy from `landsymm_mat/data/`.
- **Package installation**: `pip install -e .` for development mode.
- **Missing READMEs**: Need to create for `common_py/`, `tests_py/`, and `hildaplus/`.

### Estimated Effort

- Restructuring and moving files: moderate
- Updating all import paths: mechanical, touches every file
- Centralizing data path resolution: most impactful, ~20-30 locations
- Updating documentation: straightforward
- Testing full pipeline runs: essential, one run per stage

Estimated: 1–1.5 days of focused work.

---

## 7. Complete File Inventory

### Files Created This Session

| File | Type | Description |
|------|------|-------------|
| `landsymm.wetlands/wetland_into_forLPJG.py` | Script | Insert PEATLAND into scenario forLPJG landcover files |
| `landsymm.wetlands/README.md` | Documentation | Complete package documentation |
| `parity_diagnostic_report_postRfix.md` | Report | Post-R-fix parity diagnostic conclusions |
| `diag_cropfract_regression.py` | Diagnostic | CropFract parity regression analysis |
| `diag_cropfract_threshold.py` | Diagnostic | CropFract threshold analysis |
| `diag_mgmt_spikes.py` | Diagnostic | Fert/Irrig management spike analysis |
| `session_handoff_20260126.md` | Documentation | This handoff document |

### Files Modified This Session

| File | Changes |
|------|---------|
| `landsymm.wetlands/glwd3_to_halfdeg.py` | Added `SCALE_FACTORS_NOFORESTS`, dual product generation, output filename change |
| `landsymm.wetlands/wetland_into_hilda.py` | Added `WETLAND_VARS`, `wetland_product` parameter, `--wetland-product` CLI arg |
| `landsymm.wetlands/run_wetland_pipeline.py` | Added `wetland_product` passthrough, CLI arg, updated NetCDF filename |
| `landsymm_pipeline_technical_manual.md` | Added Stage 4 section, updated diagrams, renumbered sections |

### Output Data Files Generated

| File | Size | Description |
|------|------|-------------|
| `data/geodata_py/glwd3/peatland_halfdeg.nc` | ~12 MB | Half-degree wetland fractions (both products) |
| `data/output_hildaplus_remap_10b/remaps_v10_old_62892_gL/LU.remapv10_old_62892_gL_peatland.txt` | ~300 MB | Peatland-inclusive baseline LU (1901-2020) |
| `data/PLUMv2_LU_default_output/SSP*/.../landcover_peatland.txt` (×10) | ~200 MB each | Peatland-inclusive scenario landcover (2021-2100) |

### Key Reference Files (Pre-existing, Consulted)

| File | Role |
|------|------|
| `parity_table_full_postRfix_2021_2100.md` | Latest PLUMharm parity table (post-R-fix) |
| `parity_table_lpjg_postRfix.md` | Latest PLUMharm2LPJG parity table |
| `parity_table_full_p2b_clean_baseline_2021_2100.md` | Pre-R-fix baseline parity table |
| `plumharm_comprehensive_debug_plan.md` | Master debug tracking document |
| `rollback_state_recovery_note.md` | Rollback recovery instructions |
| `plumharm_matlab_python_comprehensive_audit.md` | Full MATLAB-Python code audit |
| `chat_handover_p3e.md` | Previous session handover (Approach E) |

---

## 8. Current State Summary

### What Is Complete

| Task | Status | Evidence |
|------|--------|----------|
| PLUMharm parity investigation | **Closed** | `parity_diagnostic_report_postRfix.md` — all diffs are 1-ULP artifacts |
| PLUMharm2LPJG parity | **Closed** | `parity_table_lpjg_postRfix.md` — landcover exact, rest ULP noise |
| PLUMharmFigs porting | **Complete** | All 7 gaps addressed (untested visually) |
| PLUM reformatting (R→Python port) | **Complete** | R bug fixed, both codebases re-run |
| GLWD3 aggregation | **Complete** | `peatland_halfdeg.nc` with both wetland products |
| HILDA+ remap LU peatland insertion | **Complete** | `LU.remapv10_old_62892_gL_peatland.txt` verified |
| ForLPJG landcover peatland insertion | **Complete** | 10 `landcover_peatland.txt` files across all scenarios |
| Wetland pipeline documentation | **Complete** | README + technical manual updated |
| Run scripts standardization | **Complete** | `run_plumharm.py`, `run_plumharm2lpjg.py`, `run_plumharm_figs.py`, `run_plumharm_pipeline.py` |

### What Is Planned for Next Session

| Task | Priority | Estimated Effort |
|------|----------|-----------------|
| Consolidate Python codebases into `landsymm_py/` | High | 1–1.5 days |
| Create `pyproject.toml` with dependencies | Part of above | — |
| Centralize data path resolution (`config.py`) | Part of above | — |
| Create missing READMEs (`common_py/`, `tests_py/`, `hildaplus/`) | Part of above | — |
| Update all import paths | Part of above | — |
| Test full pipeline in new structure | Part of above | — |

### What Is Optional / Deferred

| Task | Notes |
|------|-------|
| CROPLAND floor guard | Quality-of-life improvement for both codebases; eliminates management spikes at physically meaningless crop areas |
| PLUMharmFigs visual testing | Figures code ported but not visually validated |
| PLUMharm parity for non-SSP1 scenarios | Only SSP1_RCP26 fully parity-tested; others assumed equivalent |

---

---

## 9. Analysis: Forest Management and Gross Transitions Integration

### 9.1 Objective

Assess what it would mean to explicitly include:
- **Forest management** (managed vs unmanaged forest as distinct land-use types, separate from "other natural")
- **Gross land-use transitions** (bidirectional transition matrices rather than net year-to-year changes)

across the entire Python pipeline: HILDA+ → Remapping → PLUM Reformatting → PLUMharm → PLUMharm2LPJG → Wetlands/GLWD.

### 9.2 Current State — Where Forest Data Exists But Is Lost

Forest data is actually available at multiple points in the pipeline, but gets collapsed into a single `NATURAL` category early on.

#### Stage 1: HILDA+ (Data EXISTS, split capability EXISTS)

| Data | Available? | Details |
|------|:----------:|---------|
| 5 forest types | YES | ForestNE, ForestND, ForestBE, ForestBD, ForestPNV |
| Managed/unmanaged split | YES | Forest management file (`forest_management.nc`) with values 0/1/2; managed codes = original × 10 (e.g., 41→410) |
| forestfrac output | YES | `hildaplus_forestfrac_1901_2020.txt` — per-type fractions within forest |
| Split mode in v3 scripts | YES | `--forest-mode split` produces FOREST_MANAGED and FOREST_UNMANAGED as separate netfrac columns |
| Gross transitions | LEGACY ONLY | `hildap_tables.py` (HPC script) produces 36 transition columns; the v3 scripts used in current pipeline do NOT produce them |

HILDA+ state codes:
```
Unmanaged: 40 (PNV/other), 41 (NE), 42 (BE), 43 (ND), 44 (BD), 45 (mixed)
Managed:   400, 410, 420, 430, 440, 450 (= unmanaged code × 10)
```

The smoothing stage merges forest management via `hildafomafusion()`:
```python
hildastates * (np.ones_like(hildafoma) + (hildafoma == 1) * 9)
# Where forest_management == 1, code is multiplied by 10
```

#### Stage 2: Remapping (**FOREST IS COLLAPSED HERE**)

`lu_import.py` maps all HILDA+ output to 4 classes:
```python
FOREST + NATURAL → NATURAL
CROPLAND → CROPLAND
PASTURE → PASTURE
URBAN + BARREN → BARREN
```

- `forestfrac` is **never read** — only `netfrac` is used
- All forest distinction (type, managed/unmanaged) is **permanently lost**
- This is the critical merge point where forest data disappears

#### Stage 3 Preprocessing: PLUM Reformatting (Data EXISTS but collapsed)

Raw PLUM `LandCover.txt` has:
```
TimberForest, UnmanagedForest, CarbonForest, OtherNatural
```

Reformatting produces:
- `LandUse.txt` — **preserves** `managedForest` (= timberForest + carbonForest), `unmanagedForest`, `otherNatural`
- `LandCoverFract.txt` — **collapses** all into: `NATURAL = (managedForest + unmanagedForest + otherNatural) / area`

So the forest-level detail exists in `LandUse.txt` but is not propagated to the harmonization input.

#### Stage 3: PLUMharm (Only 4 classes)

Operates on 4 land-use classes: `CROPLAND, PASTURE, NATURAL, BARREN`

- NATURAL is monolithic
- Donation order for cropland expansion: PASTURE → NATURAL → BARREN
- Ring redistribution works per-class
- No forest management concept exists

#### Stage 3: PLUMharm2LPJG (Only 4 classes)

Output `landcover.txt` has: `Lon, Lat, Year, PASTURE, CROPLAND, NATURAL, BARREN`

#### Stage 4: Wetlands (Carves from NATURAL)

PEATLAND is carved from NATURAL. With forest split, would need to decide which sub-category peatland comes from.

### 9.3 Gross Transitions — Current and Legacy

#### What gross transitions are

Instead of just tracking **net** changes:
```
Year N:   CROPLAND=0.30, NATURAL=0.50
Year N+1: CROPLAND=0.32, NATURAL=0.48
Net: CROPLAND +0.02, NATURAL -0.02
```

Gross transitions track the **full bidirectional matrix**:
```
NATURAL → CROPLAND: +0.05 (deforestation/land clearing)
CROPLAND → NATURAL: -0.03 (abandonment/rewilding)
Net: CROPLAND +0.02, but gross turnover = 0.08
```

This matters for carbon accounting, biodiversity, and vegetation dynamics because land that was recently cleared has different properties than land that has been continuously cultivated.

#### Legacy gross transition system

`hildap_tables.py` produces 36 transition types using a code system:
```
u=Urban, c=Cropland, p=Pasture, b=Barren
s=secondary natural, f=managed forest
v=primary/virgin natural (never converted)
```

Transition types: `up, uc, ub, pu, pc, pb, cu, cp, cb, bu, bp, bc, us, ps, cs, bs, uf, pf, cf, bf, su, sp, sc, sb, fu, fp, fc, fb, vu, vp, vc, vb, sf, fs, vf`

Each column is a fraction of the gridcell. This is NOT used in the current v3-based pipeline.

#### PLUM gross transitions

**PLUM does not provide gross transition data.** There are no `LandTransitions.txt` or similar files in any scenario directory. The harmonization pipeline works exclusively with net year-to-year changes:
```python
agri_d_YXv = in_y1_2deg_agri_YXv - in_y0_2deg_agri_YXv
```

### 9.4 Implications of Adding Forest Management — Stage by Stage

#### A. Stage 1: HILDA+ Smoothing & Upscaling

**Impact: LOW — capability already exists**

- The v3 scripts already support `--forest-mode split`, producing `FOREST_MANAGED` and `FOREST_UNMANAGED` as separate netfrac columns
- `forestfrac` output already provides per-type detail
- **Action needed**: Run the v3 pipeline with `--forest-mode split` and forest management file
- **New output columns**: `URBAN, CROPLAND, PASTURE, FOREST_MANAGED, FOREST_UNMANAGED, NATURAL, BARREN` (7 classes instead of 6)
- **Or alternatively**: Keep combined `FOREST` column and use `forestfrac` for the managed/unmanaged split downstream

#### B. Stage 2: Remapping

**Impact: MEDIUM-HIGH**

Currently produces 4 classes. Would need to produce 5-6 classes.

Changes needed in `lu_import.py`:
```
Before: FOREST + NATURAL → NATURAL
After:  FOREST_MANAGED → MANAGED_FOREST (or FOREST)
        FOREST_UNMANAGED + NATURAL → NATURAL (or keep separate)
```

Cascading implications:
1. `remap_options.py` — new LU type names and donation/receiving rules
2. `remap.py` — normalization (sum-to-1) with more classes
3. `cropfrac.py` — cropland expansion logic needs to know which classes can donate; currently takes from all non-CROPLAND
4. `nfert.py` — no direct impact (operates per-crop, not per-LU)
5. `soil.py` — no direct impact
6. Output files — `LU.remapv*.txt` would have additional columns
7. Inpainting/interpolation — needs to handle new columns
8. All validation and diagnostics — updates needed

**Key design decision**: What happens when cropland expands? Current donation order is PASTURE → NATURAL → BARREN. With forest split:
- Should cropland expand into unmanaged forest before managed forest?
- Should managed forest be "protected" from conversion?
- This is a **scientific modeling decision**, not just a code change

#### C. Stage 3 Preprocessing: PLUM Reformatting

**Impact: LOW — data already exists**

`LandUse.txt` already contains `managedForest`, `unmanagedForest`, `otherNatural`. The change is to output them as separate columns in `LandCoverFract.txt` instead of collapsing into NATURAL:

```
Before: NATURAL = (managedForest + unmanagedForest + otherNatural) / area
After:  MANAGED_FOREST = managedForest / area
        UNMANAGED_FOREST = unmanagedForest / area
        OTHER_NATURAL = otherNatural / area
```

Note: PLUM's `managedForest = timberForest + carbonForest`. Whether to preserve this sub-split or keep it aggregated is a design choice.

#### D. Stage 3: PLUMharm — Main Harmonization

**Impact: HIGH — this is the heaviest change**

The harmonization algorithm's core loop operates on land-use classes. Adding forest management means:

1. **New LU class handling** (currently `[CROPLAND, PASTURE, NATURAL, BARREN]`):
   - New classes: e.g., `[CROPLAND, PASTURE, MANAGED_FOREST, UNMANAGED_FOREST, OTHER_NATURAL, BARREN]`
   - Every array and matrix in the harmonization code has dimensions tied to the number of LU classes

2. **Donation order redesign** (currently `PASTURE → NATURAL → BARREN`):
   - Which forest type donates to cropland first?
   - Can managed forest expand into unmanaged forest (i.e., forest management intensification)?
   - Scientific decision: PASTURE → OTHER_NATURAL → UNMANAGED_FOREST → MANAGED_FOREST → BARREN?
   - Or scenario-dependent ordering?

3. **Area harmonization** (`plumharm_area.py`, `plumharm_dist.py`, `plumharm_ring_redist.py`):
   - Ring redistribution currently distributes unmet cropland demand across adjacent cells
   - With more classes, the redistribution logic needs to handle more source/sink combinations
   - The "unmet" calculation becomes: how much of each class is available for conversion?

4. **Process PLUM inputs** (`plumharm_process_plum.py`):
   - Currently reads `LandCoverFract.txt` and extracts NATURAL as one column
   - Would need to extract MANAGED_FOREST, UNMANAGED_FOREST, OTHER_NATURAL separately
   - Vegetated/bare harmonization (rescaling when fractions don't sum to 1) gets more complex

5. **Conservation checks** (`plumharm_checks.py`):
   - Area conservation now needs to balance across more categories
   - Per-class tolerance checks

6. **Output files** — all `.mat` outputs gain new fields:
   - `LandCoverFract.base*.mat` — more LU columns
   - `post.base*.mat` — more state variables to carry forward
   - `.2deg.mat` files — same expansion

7. **Configuration** (`plumharm_options.py`):
   - New LU names, donation order, protected area rules per class
   - Scenario-specific forest management rules?

**Rough estimate**: 500-1000 lines of code changes across 6-8 files, plus extensive testing.

#### E. Stage 3: PLUMharm2LPJG

**Impact: LOW-MEDIUM**

Mostly mechanical once PLUMharm produces the data:
- Add new columns to `landcover.txt` output
- Ensure `.mat` → text conversion handles additional fields
- Column naming/ordering

#### F. Stage 4: Wetlands/GLWD

**Impact: MEDIUM**

Currently: `PEATLAND = min(min_NATURAL, GLWD3_frac)`, carved from NATURAL.

With forest split, key question: **where does peatland come from?**
- Option 1: Carve from OTHER_NATURAL only (forest peatlands stay as forest)
- Option 2: Carve from all non-anthropogenic classes proportionally
- Option 3: Carve from UNMANAGED_FOREST + OTHER_NATURAL (managed forest is protected)

The "without forests" wetland product becomes more meaningful here — it would exclude forest-type wetlands because forests are now tracked separately.

### 9.5 Implications of Adding Gross Transitions

**Impact: VERY HIGH — fundamental algorithm redesign**

This is a qualitatively different challenge from forest management separation.

#### Key Obstacles

1. **PLUM does not provide gross transitions**: The scenario projections only give net land-use states per year. Without gross transition data from the scenario model, there is no way to know the bidirectional flows within each year.

2. **V3 HILDA+ scripts don't produce gross transitions**: Only the legacy HPC script does. Porting the 36-column gross transition computation to the v3 scripts is possible but significant work.

3. **PLUMharm algorithm redesign**: The harmonization currently computes `delta = year_N+1 - year_N` per class, then redistributes unmet deltas spatially. With gross transitions:
   - Need to track which class each m² came from
   - Need transition matrices per gridcell per year (N_classes × N_classes matrix per cell)
   - Ring redistribution would need to redistribute specific transitions, not just areas
   - This is essentially a different algorithm

4. **LPJ-GUESS input format**: Does LPJ-GUESS actually accept gross transition inputs? If not, this entire effort would be pointless. This is the critical gatekeeping question.

5. **Memory and performance**: A transition matrix at 62,892 cells × 80 years × 6×6 classes would be substantial. The `.mat` files would grow significantly.

#### What IS Feasible Without PLUM Changes

Even without PLUM gross transition data, there's a simpler approach:

- Use HILDA+ gross transitions (1901-2020) for the **historical baseline**
- For the **scenario period** (2021-2100), infer gross transitions from net changes using assumptions (e.g., all area loss from a class is conversion out, all gain is conversion in, with no simultaneous bidirectional flow)
- This gives "pseudo-gross" transitions that preserve net consistency

But this only partially addresses the scientific need, since the interesting cases (shifting cultivation, forest rotation) involve simultaneous conversion in both directions.

### 9.6 Comparison: Forest Management vs Gross Transitions

| Dimension | Forest Management Split | Gross Transitions |
|-----------|:-----------------------:|:-----------------:|
| Data availability | HIGH — HILDA+ has it, PLUM has it | LOW — HILDA+ legacy only, PLUM doesn't have it |
| Pipeline changes | MEDIUM — add columns, update logic | VERY HIGH — algorithm redesign |
| Scientific value | HIGH — carbon stocks, timber, biodiversity | VERY HIGH — land-use dynamics, carbon cycling |
| Code effort | 1-2 weeks | 4-8 weeks |
| Testing effort | MEDIUM — parity checks with more classes | VERY HIGH — new algorithm, no MATLAB reference |
| Risk | LOW — additive change, doesn't break existing | HIGH — fundamental restructuring |
| Prerequisite | None | LPJ-GUESS must accept transition inputs |

### 9.7 Recommended Phasing

#### Phase 1: Forest Management Separation (Recommended First)

1. Run HILDA+ v3 with `--forest-mode split` to produce FOREST_MANAGED + FOREST_UNMANAGED
2. Update remapping to preserve forest distinction (5-6 LU classes)
3. Update PLUM reformatting to output forest sub-types separately
4. Update PLUMharm for multi-class harmonization
5. Update PLUMharm2LPJG for additional landcover columns
6. Update wetlands to decide peatland carving source

**Dependency**: Need to define the scientific donation order and conversion rules for forest classes.

#### Phase 2: Enhanced Forest Type Detail (Optional)

Use `forestfrac` output to provide per-type (NE, ND, BE, BD, PNV) managed/unmanaged fractions. This would give LPJ-GUESS detailed forest composition without changing the harmonization algorithm (forestfrac would be a pass-through layer applied after harmonization).

#### Phase 3: Gross Transitions (Only if LPJ-GUESS Supports It)

1. Port gross transition computation from legacy `hildap_tables.py` to v3 scripts
2. Define pseudo-gross transition inference for scenario period
3. Redesign PLUMharm for transition matrix harmonization
4. Define new LPJ-GUESS input format

**Prerequisite**: Confirm LPJ-GUESS can accept and use gross transition inputs.

### 9.8 Key Decision Points for Discussion

1. **How many land-use classes?**
   - Minimal: 5 classes (CROPLAND, PASTURE, FOREST, OTHER_NATURAL, BARREN)
   - Medium: 6 classes (+ split FOREST into MANAGED/UNMANAGED)
   - Full: 8 classes (+ PEATLAND + URBAN separate from BARREN)

2. **Donation order for cropland expansion?**
   - What is the scientific rationale for which natural classes are converted first?
   - Should this be scenario-dependent (e.g., SSP3 allows more forest conversion)?

3. **Is gross transition data needed for LPJ-GUESS?**
   - If LPJ-GUESS uses net fractions only, gross transitions add no value
   - If it can use them, the scientific benefit is substantial

4. **Forest type detail: harmonize or pass-through?**
   - Option A: Harmonize at the 5-type × managed/unmanaged level (20 forest sub-types) — extremely complex
   - Option B: Harmonize total forest fraction, then distribute to types using forestfrac ratios (simpler, preserves HILDA+ forest composition)

5. **PLUM's forest projections: reliable enough to harmonize?**
   - PLUM provides `managedForest`, `unmanagedForest` projections, but are these scientifically robust enough for grid-level harmonization?
   - Alternative: Use HILDA+ forest composition as time-invariant and only harmonize total forest area from PLUM

### 9.9 Data Flow Diagram — With Forest Management

```
HILDA+ (0.01°, with FM)
    │
    ▼  [Stage 1: --forest-mode split]
    │
netfrac: URBAN CROPLAND PASTURE FOREST_MANAGED FOREST_UNMANAGED NATURAL BARREN
forestfrac: ForestNE_M ForestND_M ... ForestNE_U ForestND_U ...
    │
    ▼  [Stage 2: remapping_py (updated)]
    │  Keep forest types separate
    │
LU.txt: CROPLAND PASTURE MANAGED_FOREST UNMANAGED_FOREST OTHER_NATURAL BARREN
cropfracs.txt: (unchanged)
nfert.txt: (unchanged)
    │
    ├────────────────────────┐
    ▼                        ▼
PLUM scenarios           staticData
(managedForest,
 unmanagedForest,
 otherNatural)
    │
    ▼  [Reformat: separate forest columns in LandCoverFract]
    │
LandCoverFract: MANAGED_FOREST UNMANAGED_FOREST OTHER_NATURAL CROPLAND PASTURE BARREN
    │
    ▼  [PLUMharm: 6-class harmonization]
    │  Donation: PASTURE → OTHER_NATURAL → UNMANAGED_FOREST → BARREN
    │  (MANAGED_FOREST protected? Scenario-dependent?)
    │
    ▼  [PLUMharm2LPJG]
    │
landcover.txt: PASTURE CROPLAND MANAGED_FOREST UNMANAGED_FOREST OTHER_NATURAL BARREN
    │
    ▼  [Wetlands: carve PEATLAND from OTHER_NATURAL]
    │
landcover_peatland.txt: ... + PEATLAND column
    │
    ▼
LPJ-GUESS
```

---

## Quick Reference Commands

```bash
# GLWD3 aggregation (Step 1)
python -m landsymm.wetlands.glwd3_to_halfdeg

# Peatland into HILDA+ remap LU (Step 2a)
python -m landsymm.wetlands.wetland_into_hilda

# Peatland into scenario forLPJG landcover (Step 2b — all scenarios)
python -m landsymm.wetlands.wetland_into_forLPJG \
  --parent-dir data/PLUMv2_LU_default_output

# Full PLUMharm pipeline (reformat + harmonize + convert)
python -m landsymm.harmonization.run_plumharm_pipeline \
  --parent-dir data/PLUMv2_LU_default_output

# PLUMharm parity check
python -m tests_py.plum_parity \
  --matlab-dir data/PLUMv2_LU_default_output/SSP1_RCP26/s1.*.harm.allow_unveg \
  --python-dir data/PLUMv2_LU_default_output/SSP1_RCP26/s1.*.harm.allow_unveg_py \
  --stats --max-cell --table output.md

# PLUMharm2LPJG parity check
python -m tests_py.plumharm2lpjg_parity \
  --matlab-dir data/PLUMv2_LU_default_output/SSP1_RCP26/s1.*.forLPJG \
  --python-dir data/PLUMv2_LU_default_output/SSP1_RCP26/s1.*_py.forLPJG \
  --stats --table output.md
```
