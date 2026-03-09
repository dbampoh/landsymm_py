"""Aggregation utilities (MATLAB-parity)."""
#
# TODO:
# - implement area-weighted aggregation (0.25->0.5, 0.5->2)
from __future__ import annotations

import numpy as np


def aggregate(array, res_in: float, res_out: float):
    """Aggregate arrays from one resolution to another (MATLAB aggregate.m)."""
    if array.ndim > 4:
        raise RuntimeError("Adapt code to work with >4 dimensions!")

    res_ratio = res_out / res_in
    if not float(res_ratio).is_integer():
        raise RuntimeError("out_res/in_res must be an integer!")
    res_ratio = int(res_ratio)

    in_array = np.asarray(array)
    tmp = np.zeros(
        (in_array.shape[0], in_array.shape[1] // res_ratio) + in_array.shape[2:],
        dtype=in_array.dtype,
    )
    for j in range(res_ratio):
        tmp += in_array[:, j::res_ratio, ...]

    out_array = np.zeros(
        (in_array.shape[0] // res_ratio, in_array.shape[1] // res_ratio) + in_array.shape[2:],
        dtype=in_array.dtype,
    )
    for i in range(res_ratio):
        out_array += tmp[i::res_ratio, ...]

    return out_array
