# `landsymm.tests` — Parity and Validation Tests

Tests for verifying MATLAB-Python parity and validating pipeline outputs.

## Running Tests

```bash
cd landsymm_py
pip install -e ".[dev]"
pytest landsymm/tests/
```

## Test Categories

- **Parity tests**: Verify that Python outputs match MATLAB outputs within floating-point tolerance
- **Regression tests**: Ensure that code changes do not alter output values
- **Validation tests**: Check output format, row sums, value ranges, etc.

## Tests in this directory

| File | What it covers |
|------|---------------|
| `test_landcover_config_parity.py` | YAML-driven HILDA+ → LPJ-GUESS mapping config (`hildaplus/config/`). Verifies that the default config and the `lpjg_v3_default` profile reproduce the historical hardcoded mapping byte-identically; that `lpjg_legacy_v1` reproduces the legacy `hildap_tables.py` policy; that `lpjg_treecrops_as_forest` loads and validates; and that the loader correctly rejects malformed configs (overlapping codes, unknown aliases, missing profiles). |
| `test_config_paths.py` | The de-hardwired path config in `landsymm/config.py`. Verifies that the `LANDSYMM_*` env vars (`DATA_DIR`, `GEODATA_DIR`, `REMAP_DIRNAME`, `REMAP_VER`, `PLUM_DIRNAME`, `MEMBER`) drive path/naming resolution; that the defaults are correct; that `harm_dirname()` and `get_remap_baseline_files()` build the expected names; and that `discover_scenarios()` finds scenario dirs (those containing a `{member}` dir). Portable: no external data needed. |
| `test_wetlands_paths.py` | `wetland_into_forLPJG` input selection (`_resolve_landcover_files`): the three modes - explicit `--landcover` files, `--flat-parent` scan (`<dir>/*/landcover.txt`), and the default forLPJG scan (`<parent>/*/<member>.*.forLPJG/landcover.txt`) - plus the `--scenarios` substring filter and missing-file error. No external data needed. |
| `test_figs.py` | Integration test for `plumharm_figs` (the diagnostic-figure module). Runs `run_plumharm_figs` on a small year window for one scenario and asserts time-series PDFs, map PNGs, and regional-statistics spreadsheets are produced. Guards the four figs-path bug fixes (None baseline-year default; `bad_base_yx` broadcast; advanced-index axis reordering at the regional-stats and management-time-series call sites). **Skipped automatically** when harmonized data is not present; point it at data via `LANDSYMM_TEST_DATA_DIR` / `LANDSYMM_TEST_REMAP_DIRNAME` / `LANDSYMM_TEST_PLUM_DIRNAME` / `LANDSYMM_TEST_SCENARIO`. |
