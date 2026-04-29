# HILDA+ → LPJ-GUESS Mapping Configuration

This directory holds the YAML-driven configuration for the
HILDA+ → LPJ-GUESS land-cover aggregation used by the
`hildap_tables_netfrac_v3*.py` upscaling scripts.

## Why a config file?

The HILDA+ → LPJ-GUESS mapping is a **policy** decision, not a fixed
property of either dataset. Different LandSyMM studies and different
LPJ-GUESS configurations have made different (defensible) choices:

- Should HILDA+ tree crops (code 23) and agroforestry (24) be lumped into
  CROPLAND, or kept distinct, or counted with FOREST?
- Should sparse/barren (code 66) be considered NATURAL or BARREN?
- Are managed (400–450) and unmanaged (40–45) forests separate LPJ-GUESS
  classes, or one combined FOREST?

Hardcoding any one answer into the upscaling scripts forces every user
to fork the codebase to change it. This config exposes the policy in a
human-readable YAML file that ships with named profiles and supports
custom overrides.

## Files

| File | Purpose |
|------|---------|
| `loader.py` | Python loader API (`load_landcover_config`, `get_categories`, `get_forest_types`, ...) |
| `lpjg_landcover_mapping.yaml` | Active default config (matches historical v3 hardcoded mapping) |
| `profiles/lpjg_v3_default.yaml` | Named version of the default (explicit `--mapping-profile lpjg_v3_default`) |
| `profiles/lpjg_legacy_v1.yaml` | Reproduces the legacy `hildap_tables.py` policy |
| `profiles/lpjg_treecrops_as_forest.yaml` | Example: tree crops & agroforestry → FOREST |

## Selecting a config

The loader's resolution order:

1. Explicit `--mapping-config /path/to/your.yaml` CLI flag (highest precedence)
2. Environment variable `LANDSYMM_LANDCOVER_CONFIG=/path/to/your.yaml`
3. `--mapping-profile <name>` CLI flag (resolves to `profiles/<name>.yaml`)
4. Environment variable `LANDSYMM_LANDCOVER_PROFILE=<name>`
5. Default: `lpjg_landcover_mapping.yaml` in this directory

Usage examples:

```bash
# Use the default (no flag needed)
python hildap_tables_netfrac_v3.py

# Use a named profile
python hildap_tables_netfrac_v3.py --mapping-profile lpjg_legacy_v1

# Use a custom config file
python hildap_tables_netfrac_v3.py --mapping-config /path/to/custom_mapping.yaml

# Via environment variable (e.g. for batch jobs)
export LANDSYMM_LANDCOVER_PROFILE=lpjg_treecrops_as_forest
python hildap_tables_netfrac_v3.py
```

## Schema

A complete config has the following top-level sections:

### `profile`
Metadata about the profile — used in error messages and logged at run start.

```yaml
profile:
  name: "my_custom_mapping"
  description: "A short, human-readable description"
  source: "HILDA+ v2"
  source_reference: "Winkler et al. (2021) doi:..."
```

### `hilda_classes`
Maps human-readable alias names to HILDA+ integer code lists. These reflect
the source data's encoding scheme. Edit only if you are pointing the
pipeline at a different source dataset (LUH2, ESA-CCI) whose integer codes
differ.

```yaml
hilda_classes:
  ocean: [0]
  urban: [11]
  annual_crops: [22]
  tree_crops: [23]
  ...
```

### `lpjguess_categories`
The genuinely policy-laden part. Each LPJ-GUESS land-cover category
references one or more `hilda_classes` aliases. The loader expands the
aliases into integer code lists at load time.

Two modes are required:

- `combined`: single FOREST column (managed + unmanaged collapsed)
- `split`: separate FOREST_MANAGED and FOREST_UNMANAGED columns

The script chooses between them at runtime based on `--forest-mode` and
data availability.

```yaml
lpjguess_categories:
  combined:
    URBAN:    [urban]
    CROPLAND: [annual_crops, tree_crops, agroforestry]
    ...
  split:
    URBAN:            [urban]
    CROPLAND:         [annual_crops, tree_crops, agroforestry]
    ...
```

**Validation**: the loader rejects a config in which the same HILDA+ code
appears in two different categories within the same mode. (E.g. assigning
`66` to both NATURAL and BARREN would be rejected — exactly one mapping
per code per mode.)

### `forest_types`
Per-PFT bucketing of HILDA+ forest codes used to write the forestfrac
output. Three variants:

- `combined` — for `forest_mode == "combined"` (each PFT bucket includes both managed and unmanaged codes)
- `managed` — for `forest_mode == "split"`, only managed codes (×10)
- `unmanaged` — for `forest_mode == "split"`, only unmanaged codes

Values are direct HILDA+ integer codes (not aliases).

```yaml
forest_types:
  combined:
    ForestNE:  [41, 410]
    ForestND:  [43, 430]
    ...
```

### `forest_management_codes`
Sentinel values used in the optional separate forest-management NetCDF
file (passed via `--fmfile`).

```yaml
forest_management_codes:
  no_forest: 0
  managed:   1
  unmanaged: 2
```

### `output`
Output formatting (rarely changed).

```yaml
output:
  dlon: 0.5
  dlat: 0.5
  precision: 7
```

## Adding a new profile

To add a new mapping profile that ships with the codebase:

1. Copy `profiles/lpjg_v3_default.yaml` to `profiles/<your_name>.yaml`.
2. Edit the `profile.name` to match the filename and update `description`.
3. Modify `lpjguess_categories` and/or `forest_types` to reflect your
   policy choice.
4. Test it loads cleanly:
   ```bash
   python -c "from hildaplus.config import load_landcover_config; \
              cfg = load_landcover_config(profile='<your_name>'); \
              print(cfg['lpjguess_categories']['combined'])"
   ```
5. Optionally add a parity test in `landsymm/tests/test_landcover_config_parity.py`.

## What this does NOT cover

- The legacy `hildap_tables.py` script (which produces gross transitions
  in addition to netfrac/forestfrac) is not yet refactored to use this
  config system. Its hardcoded policy can be inspected/copied via the
  `lpjg_legacy_v1` profile, but the script itself still uses its own
  hardcoded constants. A Phase-3 refactor would extend this config
  with a `gross_transitions` sub-schema; deferred until/unless gross
  transitions is reintegrated into the active `run_chain.sh` pipeline.
