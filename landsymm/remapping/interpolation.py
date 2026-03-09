"""Interpolation helpers (stub)."""
from __future__ import annotations

# TODO:
# - implement MATLAB inpaint_nans(method=4) equivalent


def inpaint_nans(array, method: int = 4):
    """Inpaint NaNs using MATLAB-compatible method (TODO)."""
    from landsymm.common.interpolation import inpaint_nans as _inpaint

    return _inpaint(array, method=method)
