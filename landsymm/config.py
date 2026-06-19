"""Centralized path resolution for the LandSyMM pipeline.

All data paths are resolved relative to the project root by default.
Override with environment variables to use custom data locations / naming
without editing source:

    LANDSYMM_DATA_DIR       data root (default <project_root>/data)
    LANDSYMM_GEODATA_DIR    geodata dir (default <data_dir>/geodata_py)
    LANDSYMM_REMAP_DIRNAME  remap-baseline dir name (default output_hildaplus_remap_10b)
    LANDSYMM_REMAP_VER      remap version tag (default 10_old_62892_gL)
    LANDSYMM_PLUM_DIRNAME   PLUM scenario parent dir (default PLUMv2_LU_default_output)
    LANDSYMM_MEMBER         ensemble member / median run (default s1)

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
    """Return the geodata directory.

    Checks LANDSYMM_GEODATA_DIR first (geodata often lives outside the data dir
    and varies per machine); otherwise falls back to ``<data_dir>/geodata_py``.
    """
    env_dir = os.environ.get("LANDSYMM_GEODATA_DIR")
    if env_dir:
        p = Path(env_dir)
        if not p.is_dir():
            raise FileNotFoundError(
                f"LANDSYMM_GEODATA_DIR is set to '{env_dir}' but that directory "
                "does not exist."
            )
        return p
    return get_data_dir() / "geodata_py"


def get_remap_dirname() -> str:
    """Name of the remap-baseline directory under the data dir.

    Override with LANDSYMM_REMAP_DIRNAME. Default reproduces the historical
    Stage-2 output directory name. (On this workstation the corrected baseline
    is ``output_hildaplus_remap_10b_3_py``; set the env var to use it.)
    """
    return os.environ.get("LANDSYMM_REMAP_DIRNAME") or "output_hildaplus_remap_10b"


def get_remap_ver() -> str:
    """Remap version tag used in the baseline subdir and filenames, e.g.
    ``remaps_v{ver}`` and ``LU.remapv{ver}.txt``. Override with LANDSYMM_REMAP_VER.
    """
    return os.environ.get("LANDSYMM_REMAP_VER") or "10_old_62892_gL"


def get_remap_output_dir() -> Path:
    """Return the remap output directory.

    Stage 2 (remapping) writes here; Stage 3 (harmonization) and Stage 4
    (wetlands) read the baseline from here. The directory name is configurable
    via LANDSYMM_REMAP_DIRNAME (see get_remap_dirname).
    """
    return get_data_dir() / get_remap_dirname()


def get_remap_baseline_dir() -> Path:
    """Return the versioned baseline subdir (``<remap_output_dir>/remaps_v{ver}``)."""
    return get_remap_output_dir() / f"remaps_v{get_remap_ver()}"


def get_remap_baseline_files() -> "dict[str, Path]":
    """Return the three baseline forcing files for the current remap version as a
    dict with keys ``lu``, ``cropf``, ``nfert``."""
    d = get_remap_baseline_dir()
    ver = get_remap_ver()
    return {
        "lu": d / f"LU.remapv{ver}.txt",
        "cropf": d / f"cropfracs.remapv{ver}.txt",
        "nfert": d / f"nfert.remapv{ver}.txt",
    }


def get_hildaplus_output_dir() -> Path:
    """Return the HILDA+ smoothing/upscaling output directory.

    This is where the HILDA+ pipeline (smoothing + upscaling) writes
    hildaplus_netfrac_1901_2020.txt and hildaplus_forestfrac_1901_2020.txt,
    which are consumed by the remapping step.
    """
    return get_geodata_dir() / "HILDA+" / "data" / "output"


def get_plum_dirname() -> str:
    """Name of the PLUM scenario-output (parent) directory under the data dir.

    Override with LANDSYMM_PLUM_DIRNAME (e.g. ``PLUMv2_LU_NFF_output``).
    """
    return os.environ.get("LANDSYMM_PLUM_DIRNAME") or "PLUMv2_LU_default_output"


def get_plum_output_dir() -> Path:
    """Return the PLUM scenario output (parent) directory."""
    return get_data_dir() / get_plum_dirname()


def get_member() -> str:
    """Ensemble member directory name for the forcing (the median run).

    Override with LANDSYMM_MEMBER (default ``s1``).
    """
    return os.environ.get("LANDSYMM_MEMBER") or "s1"


def harm_dirname(
    member: "str | None" = None,
    remap_ver: "str | None" = None,
    allow_unveg: bool = True,
) -> str:
    """Build the harmonized-output subdirectory name for a scenario, e.g.
    ``s1.HILDA+_remap_v10_old_62892_gL.harm.allow_unveg_py``.

    The Python pipeline uses the ``_py`` tag; PLUMharm2LPJG and PLUMharmFigs
    append ``.forLPJG`` and ``_figs`` respectively to this base name.
    """
    member = member or get_member()
    remap_ver = remap_ver or get_remap_ver()
    name = f"{member}.HILDA+_remap_v{remap_ver}.harm"
    if allow_unveg:
        name += ".allow_unveg"
    return name + "_py"


def discover_scenarios(parent_dir, member: "str | None" = None) -> "list[str]":
    """List scenario subdirectory names under ``parent_dir`` that contain a
    ``{member}`` subdir (i.e., look like PLUM scenario outputs). Sorted."""
    member = member or get_member()
    parent = Path(parent_dir)
    if not parent.is_dir():
        return []
    return sorted(
        d.name for d in parent.iterdir() if d.is_dir() and (d / member).is_dir()
    )
