"""Validation utilities (stub)."""
#
# TODO:
# - implement validation helpers for sums, NaNs, negatives
# - align tolerances with MATLAB defaults
from __future__ import annotations


def assert_sum_to_one(values, target: float = 1.0, tol: float = 1e-6):
    """Assert fractions sum to 1 within tolerance (TODO)."""
    # Pseudocode:
    # 1) Compute sums along last axis.
    # 2) Find deviations from target.
    # 3) Raise error if any deviation > tol.
    raise NotImplementedError
