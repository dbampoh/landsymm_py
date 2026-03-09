"""Gridlist utilities (stub)."""
from __future__ import annotations

# TODO:
# - implement gridlist loading with list2map/mask


def load_gridlist(path: str, **kwargs):
    """Load gridlist and list2map/mask (MATLAB-parity)."""
    from landsymm.common.gridlist import load_gridlist as _load

    return _load(path, **kwargs)
