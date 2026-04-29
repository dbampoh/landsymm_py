# scripts

This directory contains the HILDA+ smoothing and upscaling pipeline scripts,
plus the `run_chain.sh` orchestrator that ties them together.

## Active pipeline (used by `run_chain.sh`)

- `hilda_smoothing/smoothing_local.py` — Single-process smoothing of raw HILDA+ states (Gaussian filter to remove flickering artifacts).
- `hilda_smoothing/smoothing_local_parallel.py` — Multi-process version of the same.
- `hildaplus-upscaling/hildap_tables_netfrac_v3.py` — Single-process upscaling: aggregates the smoothed HILDA+ NetCDF to half-degree netfrac and forestfrac text tables. Driven by the YAML mapping config in `hildaplus/config/`.
- `hildaplus-upscaling/hildap_tables_netfrac_v3_parallel.py` — Multi-process version of the above.
- `hildaplus-upscaling/inspect_netfrac_forestfrac.py` — Post-run inspection of an already-produced netfrac text file.
- `run_chain.sh` — Pipeline orchestrator. Accepts `--mapping-config PATH` and `--mapping-profile NAME` (passed through to the upscaling scripts).

## Legacy / HPC scripts

- `hildaplus-upscaling/hildap_tables.py` — The original upscaling script. Produces gross transitions in addition to netfrac/forestfrac. Uses its own hardcoded mapping (matches the `lpjg_legacy_v1` profile in `hildaplus/config/`); not yet refactored to use the YAML loader.
- `hildaplus-upscaling/launch_tables.py` — SLURM submitter for `hildap_tables.py` on the IFU cluster (Keal).

## YAML-driven HILDA+ → LPJ-GUESS mapping

The aggregation policy used by the v3 active scripts (which HILDA+ codes go into which LPJ-GUESS class) is **configurable via YAML** rather than hardcoded. Default behavior reproduces the historical mapping. See `hildaplus/config/README.md` for the schema, available profiles, and how to add custom profiles.

## Origins

- The hildaplus-upscaling scripts originated in the group's [gitlab](https://gitlab.imk-ifu.kit.edu/lemg/belda-d/hildaplus-upscaling).
- The hilda_smoothing scripts originated in [this gitlab](https://gitlab.imk-ifu.kit.edu/lemg/wittenbrink-m/hildaplus_smoothing).   
