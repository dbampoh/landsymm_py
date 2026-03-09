"""Mapping utilities (stub)."""
#
# TODO:
# - implement vector<->map conversions with list2map
# - handle YX orientation exactly as MATLAB
from __future__ import annotations

import numpy as np


def vector_to_map(vector, shape_yx, list2map):
    """Convert vectorized cells to map (TODO)."""
    # Pseudocode:
    # 1) Allocate full map of shape_yx with NaNs.
    # 2) Assign vector values at list2map indices.
    # 3) Return map_yx.
    vector = np.asarray(vector)
    list2map = np.asarray(list2map, dtype=int)
    map_yx = np.full(shape_yx, np.nan)
    rr, cc = np.unravel_index(list2map - 1, shape_yx, order="F")
    map_yx[rr, cc] = vector
    return map_yx


def map_to_vector(map_yx, list2map):
    """Convert map to vectorized cells (TODO)."""
    # Pseudocode:
    # 1) Index map_yx by list2map.
    # 2) Return vector in list2map order.
    map_yx = np.asarray(map_yx)
    list2map = np.asarray(list2map, dtype=int)
    return map_yx.ravel(order="F")[list2map - 1]


def yxz_to_xz(in_yxz, matrix_size, list2map):
    """Convert YXz array to xz (MATLAB lpjgu_YXz_to_xz)."""
    in_yxz = np.asarray(in_yxz)
    if in_yxz.ndim != 3:
        raise RuntimeError("in_YXz must be a 3d array")
    n_cells = len(list2map)
    if matrix_size[0] != n_cells:
        raise RuntimeError(
            f"matrix_size(1) ({matrix_size[0]}) must equal length of list2map ({n_cells})"
        )
    nreps = in_yxz.shape[2]
    map_size = (in_yxz.shape[0], in_yxz.shape[1])
    list2map_all = _get_list2map_all(list2map, map_size, nreps)
    out_xz = in_yxz.ravel(order="F")[list2map_all - 1]
    out_xz = out_xz.reshape((n_cells, nreps), order="F")
    return out_xz


def xz_to_yxz(in_matrix, map_size, list2map):
    """Convert xz to YXz (MATLAB lpjgu_xz_to_YXz)."""
    in_matrix = np.asarray(in_matrix)
    if in_matrix.ndim != 2:
        raise RuntimeError("in_matrix must be a vector or 2-d array")
    if in_matrix.shape[0] != len(list2map):
        raise RuntimeError(
            f"in_matrix (size1 {in_matrix.shape[0]}) must have same first-dim size as list2map (size1 {len(list2map)})"
        )
    nreps = in_matrix.shape[1]
    list2map_all = _get_list2map_all(list2map, map_size, nreps)
    map_yxz = np.full((map_size[0], map_size[1], nreps), np.nan)
    flat = map_yxz.ravel(order="F")
    flat[list2map_all - 1] = in_matrix.ravel(order="F")
    map_yxz = flat.reshape((map_size[0], map_size[1], nreps), order="F")
    return map_yxz


def _get_list2map_all(list2map, map_size, nreps):
    n_cells = len(list2map)
    list2map_all = np.empty(n_cells * nreps, dtype=int)
    for r in range(nreps):
        i1 = r * n_cells
        iN = (r + 1) * n_cells
        list2map_all[i1:iN] = r * (map_size[0] * map_size[1]) + np.asarray(list2map)
    return list2map_all
