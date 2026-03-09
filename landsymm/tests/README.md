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
