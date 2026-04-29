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
