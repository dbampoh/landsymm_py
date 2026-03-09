"""Centralized path resolution for the LandSyMM pipeline.

All data paths are resolved relative to the project root by default.
Override with the LANDSYMM_DATA_DIR environment variable to use a
custom data location (e.g., on a shared filesystem or external drive).

Usage in any module::

    from landsymm.config import get_data_dir, get_geodata_dir

    data = get_data_dir()                    # <project_root>/data
    geo  = get_geodata_dir()                 # <project_root>/data/geodata_py
    hilda = get_geodata_dir() / "HILDA+"     # <project_root>/data/geodata_py/HILDA+
"""
from __future__ import annotations

import os
from pathlib import Path


def get_project_root() -> Path:
    """Return the landsymm_py project root directory."""
    return Path(__file__).resolve().parents[1]


def get_data_dir() -> Path:
    """Return the data directory.

    Checks LANDSYMM_DATA_DIR environment variable first; falls back to
    ``<project_root>/data/``.
    """
    env_dir = os.environ.get("LANDSYMM_DATA_DIR")
    if env_dir:
        p = Path(env_dir)
        if not p.is_dir():
            raise FileNotFoundError(
                f"LANDSYMM_DATA_DIR is set to '{env_dir}' but that directory "
                "does not exist."
            )
        return p
    return get_project_root() / "data"


def get_geodata_dir() -> Path:
    """Return the geodata directory (``data/geodata_py/``)."""
    return get_data_dir() / "geodata_py"


def get_remap_output_dir() -> Path:
    """Return the remap output directory.

    Stage 2 (remapping) writes here; Stage 3 (harmonization) and
    Stage 4 (wetlands) read baseline data from here.
    """
    return get_data_dir() / "output_hildaplus_remap_10b"


def get_hildaplus_output_dir() -> Path:
    """Return the HILDA+ smoothing/upscaling output directory.

    This is where the HILDA+ pipeline (smoothing + upscaling) writes
    hildaplus_netfrac_1901_2020.txt and hildaplus_forestfrac_1901_2020.txt,
    which are consumed by the remapping step.
    """
    return get_geodata_dir() / "HILDA+" / "data" / "output"


def get_plum_output_dir() -> Path:
    """Return the PLUM scenario output directory."""
    return get_data_dir() / "PLUMv2_LU_default_output"
