"""Loader for the HILDA+ → LPJ-GUESS land-cover mapping configuration.

This module replaces the hardcoded mapping constants that historically lived
in ``hildap_tables_netfrac_v3.py`` (and its parallel sibling and the legacy
``hildap_tables.py``). The mapping is now expressed in a YAML config file
that users can edit, swap, or override without modifying Python source.

Why this exists
---------------
The HILDA+ → LPJ-GUESS land-cover aggregation is a *policy* decision, not a
fixed property of either dataset:

- Should HILDA+ tree crops (code 23) and agroforestry (24) be lumped into
  CROPLAND, or kept distinct, or counted with FOREST?
- Should sparse/barren (code 66) be considered NATURAL or BARREN?
- Are managed (400-450) and unmanaged (40-45) forests separate LPJ-GUESS
  classes, or one combined FOREST?

Different LandSyMM studies and different LPJ-GUESS configurations have
made different (defensible) choices. Hardcoding any one choice into the
upscaling scripts forces every user to fork the codebase to change it.
This loader exposes those choices in a YAML config that ships with named
profiles (``lpjg_v3_default``, ``lpjg_legacy_v1``,
``lpjg_treecrops_as_forest``) and supports custom user-supplied configs.

Resolution order for selecting a config
---------------------------------------
1. Explicit ``config_path`` argument to ``load_landcover_config``
2. Environment variable ``LANDSYMM_LANDCOVER_CONFIG`` (full path to YAML)
3. ``profile`` argument (resolved as ``profiles/<name>.yaml``)
4. Environment variable ``LANDSYMM_LANDCOVER_PROFILE`` (profile name)
5. Default: ``lpjg_landcover_mapping.yaml`` next to this module

Companion files
---------------
- ``lpjg_landcover_mapping.yaml`` — the active default config
- ``profiles/`` — alternative named profiles (one YAML per profile)
- ``README.md`` — full schema documentation
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Optional, Union

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "PyYAML is required to load HILDA+ → LPJ-GUESS mapping configs. "
        "Install with: pip install pyyaml"
    ) from exc


CONFIG_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = CONFIG_DIR / "lpjg_landcover_mapping.yaml"
PROFILES_DIR = CONFIG_DIR / "profiles"


def _resolve_path(
    config_path: Optional[Union[str, Path]] = None,
    profile: Optional[str] = None,
) -> Path:
    """Return the YAML path to load, following the documented resolution order."""
    if config_path is not None:
        return Path(config_path)

    env_path = os.environ.get("LANDSYMM_LANDCOVER_CONFIG")
    if env_path:
        return Path(env_path)

    name = profile or os.environ.get("LANDSYMM_LANDCOVER_PROFILE")
    if name:
        candidate = PROFILES_DIR / f"{name}.yaml"
        if not candidate.exists():
            raise FileNotFoundError(
                f"Profile '{name}' not found: {candidate} does not exist. "
                f"Available profiles: {[p.stem for p in PROFILES_DIR.glob('*.yaml')]}"
            )
        return candidate

    return DEFAULT_CONFIG


def _expand_aliases(cfg: Dict) -> None:
    """Replace alias-list entries in ``lpjguess_categories`` with concrete code lists.

    In the YAML, each LPJ-GUESS category is written as a list of *aliases* into
    ``hilda_classes`` (e.g. ``[annual_crops, tree_crops]``). After this call,
    each category value is the flattened list of integer codes.
    """
    raw: Dict[str, List[int]] = cfg["hilda_classes"]
    for mode_name, categories in cfg["lpjguess_categories"].items():
        for cat_name, alias_list in list(categories.items()):
            try:
                expanded = [code for alias in alias_list for code in raw[alias]]
            except KeyError as exc:
                raise KeyError(
                    f"Profile '{cfg.get('profile', {}).get('name', '?')}', "
                    f"mode '{mode_name}', category '{cat_name}': "
                    f"unknown HILDA class alias {exc}. "
                    f"Known aliases: {sorted(raw.keys())}"
                ) from exc
            categories[cat_name] = expanded


def _validate(cfg: Dict) -> None:
    """Sanity-check the config: no overlapping categories within a mode, etc."""
    profile_name = cfg.get("profile", {}).get("name", "?")
    for mode_name, categories in cfg["lpjguess_categories"].items():
        seen: Dict[int, str] = {}
        for cat_name, codes in categories.items():
            for code in codes:
                if code in seen:
                    raise ValueError(
                        f"Profile '{profile_name}', mode '{mode_name}': "
                        f"HILDA code {code} appears in both '{seen[code]}' "
                        f"and '{cat_name}' — categories must not overlap."
                    )
                seen[code] = cat_name


def load_landcover_config(
    config_path: Optional[Union[str, Path]] = None,
    profile: Optional[str] = None,
) -> Dict:
    """Load and validate a HILDA+ → LPJ-GUESS mapping config.

    Parameters
    ----------
    config_path : str or Path, optional
        Direct path to a YAML config file. Highest precedence.
    profile : str, optional
        Named profile (resolved as ``profiles/<name>.yaml``).

    Returns
    -------
    dict
        The loaded config with ``lpjguess_categories`` already expanded from
        alias lists to integer code lists. Top-level keys include:
        ``profile``, ``hilda_classes``, ``lpjguess_categories``,
        ``forest_types``, ``forest_management_codes``, ``output``.
    """
    path = _resolve_path(config_path, profile)
    if not path.exists():
        raise FileNotFoundError(f"Land-cover config not found: {path}")

    with open(path) as f:
        cfg = yaml.safe_load(f)

    _expand_aliases(cfg)
    _validate(cfg)
    cfg.setdefault("_source_path", str(path))
    return cfg


def get_categories(
    cfg: Dict, has_fm: bool, forest_mode: str
) -> Dict[str, List[int]]:
    """Return the LPJ-GUESS land-cover categories for the given run mode.

    Equivalent of the historical ``get_lpjguess_categories(has_fm, forest_mode)``
    function, but driven by the YAML config. If ``forest_mode == "split"`` and
    forest-management info is available (``has_fm == True``), returns the
    ``split`` categories (FOREST_MANAGED + FOREST_UNMANAGED). Otherwise
    returns the ``combined`` categories (single FOREST column).
    """
    use_split = (forest_mode == "split") and has_fm
    key = "split" if use_split else "combined"
    if key not in cfg["lpjguess_categories"]:
        raise KeyError(
            f"Config has no '{key}' mode under lpjguess_categories. "
            f"Available: {list(cfg['lpjguess_categories'].keys())}"
        )
    return cfg["lpjguess_categories"][key]


def get_forest_types(cfg: Dict, variant: str = "combined") -> Dict[str, List[int]]:
    """Return the forest sub-typing dictionary used for ``forestfrac.txt``.

    Variants: ``combined`` (5 PFTs, includes both managed & unmanaged codes),
    ``managed`` (5 PFTs, only managed codes), ``unmanaged`` (5 PFTs, only
    unmanaged codes).
    """
    if variant not in cfg["forest_types"]:
        raise KeyError(
            f"Config has no '{variant}' forest_types variant. "
            f"Available: {list(cfg['forest_types'].keys())}"
        )
    return cfg["forest_types"][variant]


def get_forest_management_codes(cfg: Dict) -> Dict[str, int]:
    """Return the FM sentinel values used in a separate forest_management NetCDF."""
    return cfg.get("forest_management_codes", {"no_forest": 0, "managed": 1, "unmanaged": 2})


def get_output_settings(cfg: Dict) -> Dict:
    """Return the output formatting settings (dlon, dlat, precision)."""
    return cfg.get("output", {"dlon": 0.5, "dlat": 0.5, "precision": 7})
