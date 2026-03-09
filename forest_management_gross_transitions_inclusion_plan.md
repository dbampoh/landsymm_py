# Forest Management & Gross Land-Use Transitions — Inclusion Plan

**Date:** January 26, 2026
**Status:** Planning / Future Development
**Scope:** Full LandSyMM Python pipeline (HILDA+ → Remapping → PLUM Reformatting → PLUMharm → PLUMharm2LPJG → Wetlands/GLWD)

---

## 1. Objective

This document assesses the full implications of explicitly including:

- **Forest management** — distinguishing managed forest from unmanaged forest and from other natural land, as separate land-use types throughout the pipeline
- **Gross land-use transitions** — tracking bidirectional transition matrices (e.g., NATURAL→CROPLAND and CROPLAND→NATURAL simultaneously) rather than only net year-to-year changes

The analysis covers every pipeline stage: data availability, code changes required, scientific design decisions, effort estimates, and a recommended phased implementation.

---

## 2. Current State — Where Forest Data Exists But Gets Lost

Forest data is available at multiple points in the pipeline, but is collapsed into a single `NATURAL` category at the remapping stage (Stage 2).

### 2.1 Stage 1: HILDA+ Smoothing & Upscaling

**Forest data: FULLY AVAILABLE**

| Data | Available? | Details |
|------|:----------:|---------|
| 5 forest types | YES | ForestNE, ForestND, ForestBE, ForestBD, ForestPNV |
| Managed/unmanaged split | YES | Forest management file (`forest_management.nc`) with values 0 (no forest), 1 (managed), 2 (unmanaged) |
| forestfrac output | YES | `hildaplus_forestfrac_1901_2020.txt` — per-type fractions within total forest |
| Split mode in v3 scripts | YES | `--forest-mode split` produces FOREST_MANAGED and FOREST_UNMANAGED as separate netfrac columns |
| Gross transitions | LEGACY ONLY | `hildap_tables.py` (HPC script) produces 36 transition columns; the v3 scripts used in current pipeline do NOT produce them |

#### HILDA+ Land Cover State Codes

```
Forest (Unmanaged):  40 (PNV/other), 41 (NE), 42 (BE), 43 (ND), 44 (BD), 45 (mixed)
Forest (Managed):    400, 410, 420, 430, 440, 450 (= unmanaged code × 10)
Other land cover:    11 (Urban), 22/23/24 (Cropland), 33 (Pasture),
                     55 (Grass/shrub), 66 (Sparse/barren), 77 (Water), 99 (No data)
```

#### Forest Management Integration in Smoothing

The smoothing stage merges the forest management file into the HILDA+ states:

```python
def hildafomafusion(hildastates, hildafoma):
    """Merge HILDA+ LULC states with forest management data.
    Pixels where forest management is true are multiplied by 10."""
    return hildastates * (np.ones_like(hildafoma) + (hildafoma == 1) * 9)
```

Where `forest_management == 1`, the state code is multiplied by 10 (e.g., 41 → 410), creating integrated managed forest codes.

#### Forest Type Mapping (v3 Upscaling Scripts)

```python
FOREST_TYPES = {
    "ForestNE": [41, 410],      # Needle-leaf evergreen
    "ForestND": [43, 430],      # Deciduous needle-leaf
    "ForestBE": [42, 420],      # Broad-leaf evergreen
    "ForestBD": [44, 440],      # Broad-leaf deciduous
    "ForestPNV": [40, 45, 400, 450],  # PNV / mixed / other
}

FOREST_TYPES_MANAGED = {
    "ForestNE_MANAGED": [410], "ForestND_MANAGED": [430],
    "ForestBE_MANAGED": [420], "ForestBD_MANAGED": [440],
    "ForestPNV_MANAGED": [400, 450],
}

FOREST_TYPES_UNMANAGED = {
    "ForestNE_UNMANAGED": [41], "ForestND_UNMANAGED": [43],
    "ForestBE_UNMANAGED": [42], "ForestBD_UNMANAGED": [44],
    "ForestPNV_UNMANAGED": [40, 45],
}
```

#### HILDA+ Output Modes

| Mode | netfrac Columns | forestfrac Columns |
|------|-----------------|-------------------|
| Combined (default) | URBAN, CROPLAND, PASTURE, **FOREST**, NATURAL, BARREN | ForestNE, ForestND, ForestBE, ForestBD, ForestPNV |
| Split (`--forest-mode split`) | URBAN, CROPLAND, PASTURE, **FOREST_MANAGED, FOREST_UNMANAGED**, NATURAL, BARREN | ForestNE_MANAGED, ..., ForestPNV_UNMANAGED (10 columns) |

### 2.2 Stage 2: Remapping — **THE CRITICAL MERGE POINT**

`lu_import.py` maps all HILDA+ output to 4 classes:

```python
FOREST + NATURAL → NATURAL       # ← ALL forest distinction is lost here
CROPLAND → CROPLAND
PASTURE → PASTURE
URBAN + BARREN → BARREN
```

- `forestfrac` is **never read** — only `netfrac` is used
- All forest distinction (type, managed/unmanaged) is **permanently lost**
- This is where forest data disappears from the pipeline

### 2.3 Stage 3 Preprocessing: PLUM Reformatting

**Forest data: EXISTS but collapsed**

Raw PLUM `LandCover.txt` contains:

| Column | Description |
|--------|-------------|
| `TimberForest` | Timber production forest |
| `CarbonForest` | Carbon sequestration forest |
| `UnmanagedForest` | Unmanaged natural forest |
| `OtherNatural` | Non-forest natural land |

The reformatting script:
- Creates `managedForest = timberForest + carbonForest`
- Preserves `managedForest`, `unmanagedForest`, `otherNatural` in reformatted `LandUse.txt`
- **Collapses** all into `NATURAL = (managedForest + unmanagedForest + otherNatural) / area` in `LandCoverFract.txt`

So forest-level detail exists in `LandUse.txt` but is not propagated to the harmonization input.

### 2.4 Stage 3: PLUMharm

Operates on 4 land-use classes only: `CROPLAND, PASTURE, NATURAL, BARREN`

- NATURAL is monolithic — no forest concept
- Donation order for cropland expansion: PASTURE → NATURAL → BARREN
- Ring redistribution works per-class
- All `.mat` outputs use these 4 classes

### 2.5 Stage 3: PLUMharm2LPJG

Output `landcover.txt` has: `Lon, Lat, Year, PASTURE, CROPLAND, NATURAL, BARREN`

### 2.6 Stage 4: Wetlands/GLWD

PEATLAND is carved from NATURAL using Approach H. With forest split, would need to decide which sub-category peatland comes from.

---

## 3. Gross Transitions — Current and Legacy

### 3.1 What Gross Transitions Are

Instead of just tracking **net** changes:

```
Year N:   CROPLAND = 0.30,  NATURAL = 0.50
Year N+1: CROPLAND = 0.32,  NATURAL = 0.48
Net: CROPLAND +0.02, NATURAL −0.02
```

Gross transitions track the **full bidirectional matrix**:

```
NATURAL → CROPLAND:  +0.05  (deforestation / land clearing)
CROPLAND → NATURAL:  −0.03  (abandonment / rewilding)
Net: CROPLAND +0.02,  but gross turnover = 0.08
```

This matters for:
- **Carbon accounting**: Recently cleared land has different carbon dynamics than continuously cultivated land
- **Biodiversity**: Abandoned farmland has different ecology than never-converted natural land
- **Vegetation dynamics**: Primary vs secondary vegetation (never-converted vs re-growing)

### 3.2 Legacy Gross Transition System (HILDA+)

`hildap_tables.py` produces 36 transition types using a code system:

| Symbol | Land Type |
|--------|-----------|
| u | Urban |
| c | Cropland |
| p | Pasture |
| b | Barren |
| s | Secondary natural (previously converted, now re-growing) |
| f | Managed forest |
| v | Primary/virgin natural (never converted) |

**36 transition types**: `up, uc, ub, pu, pc, pb, cu, cp, cb, bu, bp, bc, us, ps, cs, bs, uf, pf, cf, bf, su, sp, sc, sb, fu, fp, fc, fb, vu, vp, vc, vb, sf, fs, vf`

Each column is a fraction of the gridcell transitioning between the two types per year.

**Important**: This is produced by the **legacy HPC script only** (`hildap_tables.py`). The current v3 scripts used in the pipeline do NOT produce gross transitions.

### 3.3 PLUM Gross Transitions

**PLUM does not provide gross transition data.** There are no `LandTransitions.txt` or similar files in any scenario year directory. The harmonization pipeline works exclusively with net year-to-year changes:

```python
agri_d_YXv = in_y1_2deg_agri_YXv - in_y0_2deg_agri_YXv
```

---

## 4. Forest Management Integration — Stage-by-Stage Impact

### 4.1 Stage 1: HILDA+ Smoothing & Upscaling

**Impact: LOW — capability already exists**

- The v3 scripts already support `--forest-mode split`
- `forestfrac` output already provides per-type detail
- **Action needed**: Run the v3 pipeline with `--forest-mode split` and the forest management file
- **New output columns**: `URBAN, CROPLAND, PASTURE, FOREST_MANAGED, FOREST_UNMANAGED, NATURAL, BARREN` (7 classes instead of 6)
- **Alternative**: Keep combined `FOREST` column and use `forestfrac` for the managed/unmanaged split downstream

**Changes**: Minimal — just a different invocation flag. Possibly update `run_chain.sh` defaults.

### 4.2 Stage 2: Remapping

**Impact: MEDIUM-HIGH**

Currently produces 4 classes. Would need to produce 5–6 classes.

#### Required Code Changes

| File | Change |
|------|--------|
| `lu_import.py` | Map FOREST_MANAGED → MANAGED_FOREST, FOREST_UNMANAGED + NATURAL → OTHER_NATURAL (or keep all 3 separate) |
| `remap_options.py` | New LU type names, donation rules |
| `remap.py` | Normalization (sum-to-1) with more classes |
| `cropfrac.py` | Cropland expansion logic: which classes can donate area |
| Output files | `LU.remapv*.txt` gains additional columns |
| Inpainting/interpolation | Handle new columns |
| Diagnostics/validation | Update for new class count |

`nfert.py` and `soil.py` have no direct impact (they operate per-crop or per-soil, not per-LU class).

#### Key Design Decision

What happens when cropland expands? Current donation order: PASTURE → NATURAL → BARREN.

With forest split, scientific questions arise:
- Should cropland expand into unmanaged forest before managed forest?
- Should managed forest be "protected" from conversion?
- What about OTHER_NATURAL vs forest?
- **This is a scientific modeling decision, not just a code change.**

### 4.3 Stage 3 Preprocessing: PLUM Reformatting

**Impact: LOW — data already exists**

`LandUse.txt` already contains `managedForest`, `unmanagedForest`, `otherNatural`. The change:

```
Before: NATURAL = (managedForest + unmanagedForest + otherNatural) / area
After:  MANAGED_FOREST   = managedForest   / area
        UNMANAGED_FOREST = unmanagedForest / area
        OTHER_NATURAL    = otherNatural    / area
```

Note: PLUM's `managedForest = timberForest + carbonForest`. Whether to preserve this sub-split (timber vs carbon) or keep it aggregated is a design choice.

### 4.4 Stage 3: PLUMharm — Main Harmonization

**Impact: HIGH — this is the heaviest change**

#### Changes Required

1. **New LU class handling** (currently 4 classes):
   - New set: `[CROPLAND, PASTURE, MANAGED_FOREST, UNMANAGED_FOREST, OTHER_NATURAL, BARREN]`
   - Every array and matrix in the harmonization code has dimensions tied to the number of LU classes
   - All variable naming, indexing, and slicing must be updated

2. **Donation order redesign** (currently `PASTURE → NATURAL → BARREN`):
   - Which class donates to cropland expansion first?
   - Can managed forest expand into unmanaged forest (forest management intensification)?
   - Example ordering: PASTURE → OTHER_NATURAL → UNMANAGED_FOREST → MANAGED_FOREST → BARREN
   - May need to be scenario-dependent

3. **Area harmonization** (`plumharm_area.py`, `plumharm_dist.py`, `plumharm_ring_redist.py`):
   - Ring redistribution distributes unmet cropland demand across adjacent cells
   - With more classes, the redistribution logic handles more source/sink combinations
   - The "unmet" calculation: how much of each class is available for conversion at each cell

4. **Process PLUM inputs** (`plumharm_process_plum.py`):
   - Extract MANAGED_FOREST, UNMANAGED_FOREST, OTHER_NATURAL separately from `LandCoverFract.txt`
   - Vegetated/bare harmonization (rescaling when fractions don't sum to 1) gets more complex with 6 classes

5. **Conservation checks** (`plumharm_checks.py`):
   - Area conservation must balance across more categories
   - Per-class tolerance checks

6. **Output files** — all `.mat` outputs gain new fields:
   - `LandCoverFract.base*.mat` — more LU columns
   - `post.base*.mat` — more state variables to carry forward
   - `.2deg.mat` files — same expansion

7. **Configuration** (`plumharm_options.py`):
   - New LU names, donation order, class-specific rules
   - Possibly scenario-dependent forest management policies

#### Affected Files and Modules

| Module | Change Scope |
|--------|-------------|
| `plumharm_options.py` | LU names, donation order, new parameters |
| `plumharm_import_ref.py` | Read additional LU columns from baseline |
| `plumharm_process_plum.py` | Extract forest sub-types from PLUM inputs |
| `plumharm.py` | Core loop, class indexing, output assembly |
| `plumharm_area.py` | Unmet demand calculation per class |
| `plumharm_dist.py` | Delta distribution across more classes |
| `plumharm_ring_redist.py` | Ring redistribution with more donors/receivers |
| `plumharm_mgmt.py` | Management redistribution (may need forest-specific rules) |
| `plumharm_checks.py` | Conservation checks for additional classes |

**Rough estimate**: 500–1000 lines of code changes across 8–10 files, plus extensive testing and validation.

### 4.5 Stage 3: PLUMharm2LPJG

**Impact: LOW-MEDIUM**

Mostly mechanical once PLUMharm produces the data:
- Add new columns to `landcover.txt` output
- Ensure `.mat` → text conversion handles additional fields
- Column naming/ordering

### 4.6 Stage 4: Wetlands/GLWD

**Impact: MEDIUM**

Currently: `PEATLAND = min(min_NATURAL, GLWD3_frac)`, carved from NATURAL.

With forest split, key question: **where does peatland come from?**

| Option | Source for PEATLAND | Rationale |
|--------|-------------------|-----------|
| A | OTHER_NATURAL only | Forest peatlands remain classified as forest |
| B | All non-anthropogenic proportionally | Peatland is a cross-cutting biome |
| C | UNMANAGED_FOREST + OTHER_NATURAL | Managed forest is protected from conversion |

The "without forests" wetland product becomes more meaningful: it would exclude forest-type wetlands because forests are now tracked separately as a land-use class.

---

## 5. Gross Transitions Integration — Impact Assessment

**Impact: VERY HIGH — fundamental algorithm redesign**

This is a qualitatively different challenge from forest management separation.

### 5.1 Key Obstacles

1. **PLUM does not provide gross transitions**: Scenario projections give only net land-use states per year. Without gross transition data from the scenario model, there is no way to know bidirectional flows within each year.

2. **V3 HILDA+ scripts don't produce gross transitions**: Only the legacy HPC script (`hildap_tables.py`) does. Porting the 36-column gross transition computation to the v3 scripts is possible but significant work.

3. **PLUMharm algorithm redesign**: The harmonization currently computes `delta = year_N+1 − year_N` per class, then redistributes unmet deltas spatially. With gross transitions:
   - Need to track which class each m² came from
   - Need transition matrices per gridcell per year (N_classes × N_classes)
   - Ring redistribution must redistribute specific transitions, not just areas
   - This is essentially a **different algorithm**

4. **LPJ-GUESS input format**: Does LPJ-GUESS accept gross transition inputs? If not, this effort has no downstream consumer. **This is the critical gatekeeping question.**

5. **Memory and performance**: A transition matrix at 62,892 cells × 80 years × 6×6 classes would be substantial. The `.mat` files would grow significantly.

### 5.2 What IS Feasible Without PLUM Changes

Even without PLUM gross transition data, a simpler approach exists:

- Use HILDA+ gross transitions (1901–2020) for the **historical baseline**
- For the **scenario period** (2021–2100), infer gross transitions from net changes using assumptions:
  - All area loss from a class = conversion out
  - All area gain = conversion in
  - No simultaneous bidirectional flow assumed
- This gives "pseudo-gross" transitions that preserve net consistency

**Limitation**: This only partially addresses the scientific need. The interesting cases (shifting cultivation, forest rotation) involve simultaneous conversion in both directions, which cannot be inferred from net changes alone.

### 5.3 Implementation Steps (If Pursued)

1. Port gross transition computation from legacy `hildap_tables.py` to v3 scripts
2. Update remapping to propagate transition matrices through the baseline
3. Define pseudo-gross transition inference rules for the scenario period
4. Redesign PLUMharm for transition matrix harmonization
5. Define new LPJ-GUESS input format for gross transitions
6. Extensive testing (no MATLAB reference exists for this)

---

## 6. Comparison: Forest Management vs Gross Transitions

| Dimension | Forest Management Split | Gross Transitions |
|-----------|:-----------------------:|:-----------------:|
| **Data availability** | HIGH — HILDA+ has it, PLUM has it | LOW — HILDA+ legacy only, PLUM doesn't provide it |
| **Pipeline changes** | MEDIUM — add columns, update logic | VERY HIGH — algorithm redesign |
| **Scientific value** | HIGH — carbon stocks, timber, biodiversity | VERY HIGH — land-use dynamics, carbon cycling |
| **Code effort** | 1–2 weeks | 4–8 weeks |
| **Testing effort** | MEDIUM — parity checks with more classes | VERY HIGH — new algorithm, no MATLAB reference |
| **Risk** | LOW — additive change, doesn't break existing | HIGH — fundamental restructuring |
| **Prerequisite** | None | LPJ-GUESS must accept transition inputs |

---

## 7. Recommended Phasing

### Phase 1: Forest Management Separation (Recommended First)

**Estimated effort: 1–2 weeks**

1. Run HILDA+ v3 with `--forest-mode split` to produce FOREST_MANAGED + FOREST_UNMANAGED
2. Update remapping (`lu_import.py`, `remap.py`, `remap_options.py`) to preserve forest distinction (5–6 LU classes)
3. Update PLUM reformatting (`reformat_plum_gridded.py`) to output forest sub-types separately in `LandCoverFract.txt`
4. Update PLUMharm for multi-class harmonization (donation order, ring redistribution, conservation checks)
5. Update PLUMharm2LPJG for additional landcover columns
6. Update wetlands to decide peatland carving source
7. End-to-end testing

**Prerequisite**: Define the scientific donation order and conversion rules for forest classes.

### Phase 2: Enhanced Forest Type Detail (Optional)

**Estimated effort: 3–5 days**

Use `forestfrac` output to provide per-type (NE, ND, BE, BD, PNV) managed/unmanaged fractions. This would give LPJ-GUESS detailed forest composition without changing the harmonization algorithm — forestfrac would be a pass-through layer applied after harmonization.

Two approaches:
- **Option A**: Harmonize at the 5-type × managed/unmanaged level (10 forest sub-types) — extremely complex
- **Option B**: Harmonize total managed/unmanaged forest fraction, then distribute to types using HILDA+ forestfrac ratios (simpler, preserves HILDA+ forest composition)

Option B is strongly recommended.

### Phase 3: Gross Transitions (Only if LPJ-GUESS Supports It)

**Estimated effort: 4–8 weeks**

1. Port gross transition computation from legacy `hildap_tables.py` to v3 scripts
2. Define pseudo-gross transition inference for scenario period (2021–2100)
3. Redesign PLUMharm for transition matrix harmonization
4. Define new LPJ-GUESS input format

**Prerequisite**: Confirm LPJ-GUESS can accept and use gross transition inputs. If LPJ-GUESS only uses net land-use fractions, gross transitions add no value to the pipeline.

---

## 8. Key Decision Points

These scientific and design decisions must be resolved before implementation:

### 8.1 Number of Land-Use Classes

| Level | Classes | Notes |
|-------|---------|-------|
| Minimal (current + forest) | 5: CROPLAND, PASTURE, FOREST, OTHER_NATURAL, BARREN | Simplest upgrade; forest is one block |
| Medium (managed/unmanaged) | 6: CROPLAND, PASTURE, MANAGED_FOREST, UNMANAGED_FOREST, OTHER_NATURAL, BARREN | Distinguishes forest management intensity |
| Full | 8: + PEATLAND + URBAN (separate from BARREN) | Maximum granularity |

### 8.2 Donation Order for Cropland Expansion

When cropland area increases, which classes lose area first?

Current: `PASTURE → NATURAL → BARREN`

With forest split, possible orderings:
- `PASTURE → OTHER_NATURAL → UNMANAGED_FOREST → MANAGED_FOREST → BARREN` (managed forest last — protected)
- `PASTURE → OTHER_NATURAL → UNMANAGED_FOREST → BARREN` (managed forest never donates)
- Scenario-dependent: SSP3 (regional rivalry) may allow more forest conversion than SSP1 (sustainability)

**This is a scientific modeling decision.**

### 8.3 LPJ-GUESS Compatibility

- Does LPJ-GUESS accept > 4 land-use classes in its landcover input?
- Does LPJ-GUESS accept gross transition inputs?
- What column names / format does LPJ-GUESS expect?

**These must be confirmed before any implementation.**

### 8.4 Forest Type Detail Strategy

| Strategy | Approach | Complexity |
|----------|----------|:----------:|
| Harmonize per-type | Treat each of 10 forest sub-types as a separate LU class | VERY HIGH |
| Pass-through forestfrac | Harmonize total forest fraction; distribute to types using HILDA+ forestfrac ratios | LOW |
| PLUM-driven split | Use PLUM's managedForest/unmanagedForest projections directly | MEDIUM |

Recommendation: **Pass-through forestfrac** (harmonize total managed and unmanaged forest areas; distribute to NE/ND/BE/BD/PNV using the HILDA+ baseline forestfrac composition, held time-invariant or evolving with a simple model).

### 8.5 PLUM Forest Projections Quality

PLUM provides `managedForest`, `unmanagedForest` projections per scenario. Questions:
- Are these scientifically robust enough for grid-level (0.5°) harmonization?
- Do they align with HILDA+ definitions of managed/unmanaged?
- Alternative: Use HILDA+ forest composition as time-invariant and only harmonize **total forest area** from PLUM

---

## 9. Data Flow Diagrams

### 9.1 Current Pipeline (Net, 4 classes)

```
HILDA+ (0.01°) → netfrac: URBAN CROPLAND PASTURE FOREST NATURAL BARREN
    │
    ▼  [Remapping: FOREST+NATURAL → NATURAL, URBAN → BARREN]
    │
LU.txt: CROPLAND PASTURE NATURAL BARREN  (4 classes)
    │
    ▼  [PLUMharm: net harmonization on 4 classes]
    │
landcover.txt: PASTURE CROPLAND NATURAL BARREN
    │
    ▼  [Wetlands: PEATLAND from NATURAL]
    │
landcover_peatland.txt: ... + PEATLAND
```

### 9.2 With Forest Management (Net, 6 classes)

```
HILDA+ (0.01°, with FM, --forest-mode split)
    │
    ▼
netfrac: URBAN CROPLAND PASTURE FOREST_MANAGED FOREST_UNMANAGED NATURAL BARREN
forestfrac: ForestNE_M ForestND_M ... ForestNE_U ForestND_U ...
    │
    ▼  [Remapping: keep forest types separate]
    │
LU.txt: CROPLAND PASTURE MANAGED_FOREST UNMANAGED_FOREST OTHER_NATURAL BARREN  (6 classes)
cropfracs.txt: (unchanged)
nfert.txt: (unchanged)
    │
    ├────────────────────────────────────────┐
    ▼                                        ▼
PLUM scenarios                           staticData
(managedForest, unmanagedForest,
 otherNatural)
    │
    ▼  [Reformat: separate forest columns]
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

### 9.3 With Gross Transitions (Full, 6 classes)

```
HILDA+ (0.01°, with FM)
    │
    ▼  [v3 scripts + gross transition port]
    │
netfrac: 6-class net fractions (as above)
grossfrac: 36-column transition matrix per gridcell per year
    │
    ▼  [Remapping: propagate transitions]
    │
LU.txt: 6-class net fractions
LU_transitions.txt: per-class transition matrices
    │
    ├───────────────────────────────────┐
    ▼                                   ▼
PLUM scenarios                      staticData
(net states only — no gross)
    │
    ▼  [Reformat + pseudo-gross inference]
    │
LandCoverFract: 6-class fractions
Transitions: inferred from net differences (pseudo-gross)
    │
    ▼  [PLUMharm: transition-matrix harmonization]
    │  Redistribute specific transitions, not just areas
    │
    ▼  [PLUMharm2LPJG]
    │
landcover.txt: 6-class fractions
transitions.txt: per-class transition fractions
    │
    ▼
LPJ-GUESS (must accept transition input format)
```

---

## 10. Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|:----------:|:------:|------------|
| LPJ-GUESS doesn't accept forest split | LOW | HIGH | Confirm before starting Phase 1 |
| LPJ-GUESS doesn't accept gross transitions | MEDIUM | VERY HIGH | Confirm before starting Phase 3 |
| PLUM forest projections don't align with HILDA+ | MEDIUM | MEDIUM | Use HILDA+ composition as fallback |
| Donation order choice impacts results significantly | HIGH | MEDIUM | Sensitivity analysis across orderings |
| Testing complexity with more classes | HIGH | MEDIUM | Phased testing, parity checks per stage |
| Performance degradation with transition matrices | MEDIUM | LOW | Profile and optimize as needed |

---

## 11. Summary

Forest management separation is a tractable, valuable, and well-supported enhancement. The data exists at both ends of the pipeline (HILDA+ and PLUM), the HILDA+ v3 scripts already support split mode, and the changes cascade predictably through the pipeline.

Gross transitions are a fundamentally different and much larger undertaking. The critical blocker is that PLUM does not provide gross transition data, and the prerequisite of LPJ-GUESS support must be confirmed first. If pursued, it requires an algorithm redesign in PLUMharm and a new output format.

The recommended path is: **Phase 1 first** (forest management, 1–2 weeks), then evaluate whether Phase 2 (forest type detail) or Phase 3 (gross transitions) provides more scientific value for the modeling objectives.

---

## References

- Lehner, B. & Döll, P. (2004). Development and validation of a global database of lakes, reservoirs and wetlands. *Journal of Hydrology*, 296, 1-22.
- HILDA+ v2.0 land-use/land-cover classification system
- LPJ-GUESS model documentation (for input format confirmation)
