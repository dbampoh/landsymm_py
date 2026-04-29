"""Parity tests for the YAML-driven HILDA+ → LPJ-GUESS land-cover mapping.

These tests guarantee that the YAML default config (and the named
``lpjg_v3_default`` profile) reproduce the historical hardcoded mapping
that lived in ``hildap_tables_netfrac_v3.py`` exactly. If anything in the
default YAML drifts from the original constants, these tests will catch it.

A separate test verifies that ``lpjg_legacy_v1`` reproduces the historical
``hildap_tables.py`` mapping (different policy from v3).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Make `hildaplus.config` importable (project root is one level above this
# tests/ directory's parent, which is the `landsymm/` package root).
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from hildaplus.config import (  # noqa: E402  (import after sys.path tweak)
    load_landcover_config,
    get_categories,
    get_forest_types,
    get_forest_management_codes,
    get_output_settings,
)


# ---------------------------------------------------------------------------
# Reference values copied verbatim from the historical hardcoded constants
# in hildap_tables_netfrac_v3.py
# ---------------------------------------------------------------------------

V3_LPJG_CATEGORIES_COMBINED = {
    "URBAN":    [11],
    "CROPLAND": [22, 23, 24],
    "PASTURE":  [33],
    "FOREST":   [40, 41, 42, 43, 44, 45, 400, 410, 420, 430, 440, 450],
    "NATURAL":  [55, 66],
    "BARREN":   [0, 77, 99],
}

V3_LPJG_CATEGORIES_SPLIT = {
    "URBAN":            [11],
    "CROPLAND":         [22, 23, 24],
    "PASTURE":          [33],
    "FOREST_MANAGED":   [400, 410, 420, 430, 440, 450],
    "FOREST_UNMANAGED": [40, 41, 42, 43, 44, 45],
    "NATURAL":          [55, 66],
    "BARREN":           [0, 77, 99],
}

V3_FOREST_TYPES_COMBINED = {
    "ForestNE":  [41, 410],
    "ForestND":  [43, 430],
    "ForestBE":  [42, 420],
    "ForestBD":  [44, 440],
    "ForestPNV": [40, 45, 400, 450],
}

V3_FOREST_TYPES_MANAGED = {
    "ForestNE_MANAGED":  [410],
    "ForestND_MANAGED":  [430],
    "ForestBE_MANAGED":  [420],
    "ForestBD_MANAGED":  [440],
    "ForestPNV_MANAGED": [400, 450],
}

V3_FOREST_TYPES_UNMANAGED = {
    "ForestNE_UNMANAGED":  [41],
    "ForestND_UNMANAGED":  [43],
    "ForestBE_UNMANAGED":  [42],
    "ForestBD_UNMANAGED":  [44],
    "ForestPNV_UNMANAGED": [40, 45],
}

# Reference values for the legacy hildap_tables.py policy.
# Differences from v3:
#   - sparse_barren (66) goes to BARREN, not NATURAL
#   - CROPLAND has only [22] (no tree_crops 23 or agroforestry 24)
LEGACY_LPJG_CATEGORIES_COMBINED = {
    "URBAN":    [11],
    "CROPLAND": [22],
    "PASTURE":  [33],
    "FOREST":   [40, 41, 42, 43, 44, 45, 400, 410, 420, 430, 440, 450],
    "NATURAL":  [55],
    "BARREN":   [0, 77, 99, 66],
}


# ---------------------------------------------------------------------------
# Parity tests for the default config (== lpjg_v3_default profile)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def default_cfg():
    return load_landcover_config()


def test_default_profile_name(default_cfg):
    assert default_cfg["profile"]["name"] == "lpjg_v3_default"


@pytest.mark.parametrize(
    "forest_mode,has_fm,expected",
    [
        ("combined", False, V3_LPJG_CATEGORIES_COMBINED),
        ("combined", True,  V3_LPJG_CATEGORIES_COMBINED),
        ("split",    True,  V3_LPJG_CATEGORIES_SPLIT),
        ("split",    False, V3_LPJG_CATEGORIES_COMBINED),  # split falls back to combined when no FM
    ],
)
def test_default_categories_match_v3_constants(default_cfg, forest_mode, has_fm, expected):
    actual = get_categories(default_cfg, has_fm=has_fm, forest_mode=forest_mode)
    assert actual == expected


@pytest.mark.parametrize(
    "variant,expected",
    [
        ("combined",  V3_FOREST_TYPES_COMBINED),
        ("managed",   V3_FOREST_TYPES_MANAGED),
        ("unmanaged", V3_FOREST_TYPES_UNMANAGED),
    ],
)
def test_default_forest_types_match_v3_constants(default_cfg, variant, expected):
    actual = get_forest_types(default_cfg, variant)
    assert actual == expected


def test_default_fm_codes(default_cfg):
    assert get_forest_management_codes(default_cfg) == {
        "no_forest": 0,
        "managed":   1,
        "unmanaged": 2,
    }


def test_default_output_settings(default_cfg):
    out = get_output_settings(default_cfg)
    assert out["dlon"] == 0.5
    assert out["dlat"] == 0.5
    assert out["precision"] == 7


# ---------------------------------------------------------------------------
# Named profile tests
# ---------------------------------------------------------------------------

def test_v3_default_profile_matches_default_config():
    """The 'lpjg_v3_default' named profile should be byte-equivalent to the default."""
    via_profile = load_landcover_config(profile="lpjg_v3_default")
    via_default = load_landcover_config()
    # Drop ``_source_path`` (differs by definition) before comparing.
    via_profile.pop("_source_path", None)
    via_default.pop("_source_path", None)
    assert via_profile == via_default


def test_legacy_profile_matches_legacy_constants():
    cfg = load_landcover_config(profile="lpjg_legacy_v1")
    actual = get_categories(cfg, has_fm=False, forest_mode="combined")
    assert actual == LEGACY_LPJG_CATEGORIES_COMBINED


def test_treecrops_as_forest_profile_loads_and_validates():
    """Smoke test: alternative profile should load and validate cleanly."""
    cfg = load_landcover_config(profile="lpjg_treecrops_as_forest")
    cats = get_categories(cfg, has_fm=False, forest_mode="combined")
    # Tree crops (23) and agroforestry (24) should now be in FOREST, not CROPLAND
    assert 23 in cats["FOREST"]
    assert 24 in cats["FOREST"]
    assert 23 not in cats["CROPLAND"]
    assert 24 not in cats["CROPLAND"]
    assert cats["CROPLAND"] == [22]


# ---------------------------------------------------------------------------
# Loader behavior tests
# ---------------------------------------------------------------------------

def test_unknown_profile_raises():
    with pytest.raises(FileNotFoundError, match="Profile 'no_such_profile' not found"):
        load_landcover_config(profile="no_such_profile")


def test_overlap_validation(tmp_path):
    """A config that puts the same code in two categories should be rejected."""
    bad_yaml = tmp_path / "bad.yaml"
    bad_yaml.write_text(
        """
profile:
  name: bad_overlap
hilda_classes:
  ocean: [0]
  urban: [11]
  pasture: [33]
  forest_unmanaged: [40, 41]
  forest_managed: [400, 410]
  grassland_shrubland: [55]
  sparse_barren: [66]
  water: [77]
  nodata: [99]
  annual_crops: [22]
lpjguess_categories:
  combined:
    URBAN:    [urban]
    CROPLAND: [annual_crops, urban]    # 11 appears in both URBAN and CROPLAND
    PASTURE:  [pasture]
    FOREST:   [forest_unmanaged, forest_managed]
    NATURAL:  [grassland_shrubland, sparse_barren]
    BARREN:   [ocean, water, nodata]
forest_types:
  combined:
    ForestNE: [41, 410]
forest_management_codes:
  no_forest: 0
  managed: 1
  unmanaged: 2
output:
  dlon: 0.5
  dlat: 0.5
  precision: 7
"""
    )
    with pytest.raises(ValueError, match="HILDA code 11 appears in both"):
        load_landcover_config(config_path=bad_yaml)


def test_unknown_alias_raises(tmp_path):
    """Referencing a non-existent hilda_classes alias should fail loudly."""
    bad_yaml = tmp_path / "bad_alias.yaml"
    bad_yaml.write_text(
        """
profile:
  name: bad_alias
hilda_classes:
  ocean: [0]
lpjguess_categories:
  combined:
    BARREN: [does_not_exist]
forest_types:
  combined:
    ForestNE: [41]
forest_management_codes:
  no_forest: 0
  managed: 1
  unmanaged: 2
output:
  dlon: 0.5
  dlat: 0.5
  precision: 7
"""
    )
    with pytest.raises(KeyError, match="unknown HILDA class alias"):
        load_landcover_config(config_path=bad_yaml)
