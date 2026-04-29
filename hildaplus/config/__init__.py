"""Configuration package for the HILDA+ → LPJ-GUESS upscaling pipeline.

This package provides a YAML-driven mapping configuration that lets users
change the HILDA+ → LPJ-GUESS land-cover aggregation policy without
editing Python source code. See `loader.py` for the API and
`README.md` for the schema.
"""
from .loader import (
    load_landcover_config,
    get_categories,
    get_forest_types,
    get_forest_management_codes,
    get_output_settings,
    DEFAULT_CONFIG,
    PROFILES_DIR,
)

__all__ = [
    "load_landcover_config",
    "get_categories",
    "get_forest_types",
    "get_forest_management_codes",
    "get_output_settings",
    "DEFAULT_CONFIG",
    "PROFILES_DIR",
]
