"""Interpolation utilities."""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsqr, spsolve


def inpaint_nans(*args, **kwargs):
    """Inpaint NaNs using MATLAB-compatible method."""
    array = args[0]
    method = kwargs.get("method", 0) if len(args) < 2 else args[1]

    a = np.asarray(array, dtype=float)
    n, m = a.shape
    a_flat = a.ravel(order="F")
    nan_mask = np.isnan(a_flat)
    if not np.any(nan_mask):
        return a.copy()

    # 0-based indices for Python
    nan_list = np.where(nan_mask)[0]
    known_list = np.where(~nan_mask)[0]
    nan_count = len(nan_list)

    # convert to (row,col) 0-based
    rows, cols = np.unravel_index(nan_list, (n, m), order="F")
    nan_list_rc = np.column_stack([nan_list, rows, cols])

    if method not in (0, 1, 2, 3, 4, 5):
        raise ValueError("If supplied, method must be one of: {0,1,2,3,4,5}.")

    nm = n * m

    if method == 0:
        if (m == 1) or (n == 1):
            work_list = nan_list.copy()
            work_list = np.unique(np.concatenate([work_list, work_list - 1, work_list + 1]))
            work_list = work_list[(work_list > 0) & (work_list < (nm - 1))]
            nw = work_list.size
            if nw == 0:
                return a.copy()
            u = np.arange(nw)
            cols_idx = np.vstack([work_list - 1, work_list, work_list + 1]).T
            rows_idx = np.repeat(u, 3)
            data = np.tile([1.0, -2.0, 1.0], nw)
            fda = sparse.coo_matrix((data, (rows_idx, cols_idx.ravel())), shape=(nw, nm)).tocsr()
        else:
            talks_to = np.array([[-1, 0], [0, -1], [1, 0], [0, 1]])
            neighbors_list = _identify_neighbors(n, m, nan_list_rc, talks_to)
            all_list = np.vstack([nan_list_rc, neighbors_list])

            L = (all_list[:, 1] > 0) & (all_list[:, 1] < (n - 1))
            fda = sparse.csr_matrix((nm, nm))
            if np.any(L):
                idx = all_list[L, 0]
                rows_idx = np.repeat(idx, 3)
                cols_idx = np.vstack([idx - 1, idx, idx + 1]).T.ravel()
                data = np.tile([1.0, -2.0, 1.0], idx.size)
                fda = sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()
            else:
                fda = sparse.csr_matrix((nm, nm))

            L = (all_list[:, 2] > 0) & (all_list[:, 2] < (m - 1))
            if np.any(L):
                idx = all_list[L, 0]
                rows_idx = np.repeat(idx, 3)
                cols_idx = np.vstack([idx - n, idx, idx + n]).T.ravel()
                data = np.tile([1.0, -2.0, 1.0], idx.size)
                fda = fda + sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()

        rhs = -fda[:, known_list].dot(a_flat[known_list])
        k = np.where(fda[:, nan_list].getnnz(axis=1) > 0)[0]
        A_nan = fda[k][:, nan_list]
        sol = lsqr(A_nan, rhs[k])[0]
        b_flat = a_flat.copy()
        b_flat[nan_list] = sol
        return b_flat.reshape((n, m), order="F")

    if method == 1:
        if (m == 1) or (n == 1):
            u = np.arange(nm - 2)
            rows_idx = np.repeat(u, 3)
            cols_idx = np.vstack([u, u + 1, u + 2]).T.ravel()
            data = np.tile([1.0, -2.0, 1.0], nm - 2)
            fda = sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm - 2, nm)).tocsr()
        else:
            i, j = np.mgrid[1 : n - 1, 0:m]
            ind = (i + j * n).ravel()
            fda = sparse.coo_matrix(
                (
                    np.tile([1.0, -2.0, 1.0], ind.size),
                    (np.repeat(ind, 3), np.vstack([ind - 1, ind, ind + 1]).T.ravel()),
                ),
                shape=(nm, nm),
            ).tocsr()

            i, j = np.mgrid[0:n, 1 : m - 1]
            ind = (i + j * n).ravel()
            fda = fda + sparse.coo_matrix(
                (
                    np.tile([1.0, -2.0, 1.0], ind.size),
                    (np.repeat(ind, 3), np.vstack([ind - n, ind, ind + n]).T.ravel()),
                ),
                shape=(nm, nm),
            ).tocsr()

        rhs = -fda[:, known_list].dot(a_flat[known_list])
        k = np.where(fda[:, nan_list].getnnz(axis=1) > 0)[0]
        A_nan = fda[k][:, nan_list]
        sol = lsqr(A_nan, rhs[k])[0]
        b_flat = a_flat.copy()
        b_flat[nan_list] = sol
        return b_flat.reshape((n, m), order="F")

    if method == 2:
        if (m == 1) or (n == 1):
            raise ValueError("Method 2 has problems for vector input. Please use another method.")

        fda = sparse.csr_matrix((nm, nm))
        L = (nan_list_rc[:, 1] > 0) & (nan_list_rc[:, 1] < (n - 1))
        if np.any(L):
            idx = nan_list_rc[L, 0]
            rows_idx = np.repeat(idx, 3)
            cols_idx = np.vstack([idx - 1, idx, idx + 1]).T.ravel()
            data = np.tile([1.0, -2.0, 1.0], idx.size)
            fda = sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()
        else:
            fda = sparse.csr_matrix((nm, nm))

        L = (nan_list_rc[:, 2] > 0) & (nan_list_rc[:, 2] < (m - 1))
        if np.any(L):
            idx = nan_list_rc[L, 0]
            rows_idx = np.repeat(idx, 3)
            cols_idx = np.vstack([idx - n, idx, idx + n]).T.ravel()
            data = np.tile([1.0, -2.0, 1.0], idx.size)
            fda = fda + sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()

        if 0 in nan_list:
            fda[0, [0, 1, n]] = [-2.0, 1.0, 1.0]
        if (n - 1) in nan_list:
            fda[n - 1, [n - 1, n - 2, n + (n - 1)]] = [-2.0, 1.0, 1.0]
        if (nm - n) in nan_list:
            fda[nm - n, [nm - n, nm - n + 1, nm - (n + 1)]] = [-2.0, 1.0, 1.0]
        if (nm - 1) in nan_list:
            fda[nm - 1, [nm - 1, nm - 2, nm - (n + 1)]] = [-2.0, 1.0, 1.0]

        rhs = -fda[:, known_list].dot(a_flat[known_list])
        b_flat = a_flat.copy()
        b_flat[nan_list] = spsolve(fda[nan_list][:, nan_list], rhs[nan_list])
        return b_flat.reshape((n, m), order="F")

    if method == 3:
        talks_to = np.array(
            [
                [-2, 0],
                [-1, -1],
                [-1, 0],
                [-1, 1],
                [0, -2],
                [0, -1],
                [0, 1],
                [0, 2],
                [1, -1],
                [1, 0],
                [1, 1],
                [2, 0],
            ]
        )
        neighbors_list = _identify_neighbors(n, m, nan_list_rc, talks_to)
        all_list = np.vstack([nan_list_rc, neighbors_list])

        fda = sparse.csr_matrix((nm, nm))
        L = (
            (all_list[:, 1] >= 2)
            & (all_list[:, 1] <= (n - 3))
            & (all_list[:, 2] >= 2)
            & (all_list[:, 2] <= (m - 3))
        )
        if np.any(L):
            idx = all_list[L, 0]
            rows_idx = np.repeat(idx, 13)
            cols_offsets = np.array(
                [-2 * n, -n - 1, -n, -n + 1, -2, -1, 0, 1, 2, n - 1, n, n + 1, 2 * n]
            )
            cols_idx = (idx[:, None] + cols_offsets[None, :]).ravel()
            data = np.tile([1.0, 2.0, -8.0, 2.0, 1.0, -8.0, 20.0, -8.0, 1.0, 2.0, -8.0, 2.0, 1.0], idx.size)
            fda = sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()

        L = (
            ((all_list[:, 1] == 1) | (all_list[:, 1] == (n - 2)))
            & (all_list[:, 2] >= 1)
            & (all_list[:, 2] <= (m - 2))
        ) | (
            ((all_list[:, 2] == 1) | (all_list[:, 2] == (m - 2)))
            & (all_list[:, 1] >= 1)
            & (all_list[:, 1] <= (n - 2))
        )
        if np.any(L):
            idx = all_list[L, 0]
            rows_idx = np.repeat(idx, 5)
            cols_offsets = np.array([-n, -1, 0, 1, n])
            cols_idx = (idx[:, None] + cols_offsets[None, :]).ravel()
            data = np.tile([1.0, 1.0, -4.0, 1.0, 1.0], idx.size)
            fda = fda + sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()

        L = ((all_list[:, 1] == 0) | (all_list[:, 1] == (n - 1))) & (
            (all_list[:, 2] >= 1) & (all_list[:, 2] <= (m - 2))
        )
        if np.any(L):
            idx = all_list[L, 0]
            rows_idx = np.repeat(idx, 3)
            cols_offsets = np.array([-n, 0, n])
            cols_idx = (idx[:, None] + cols_offsets[None, :]).ravel()
            data = np.tile([1.0, -2.0, 1.0], idx.size)
            fda = fda + sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()

        L = ((all_list[:, 2] == 0) | (all_list[:, 2] == (m - 1))) & (
            (all_list[:, 1] >= 1) & (all_list[:, 1] <= (n - 2))
        )
        if np.any(L):
            idx = all_list[L, 0]
            rows_idx = np.repeat(idx, 3)
            cols_offsets = np.array([-1, 0, 1])
            cols_idx = (idx[:, None] + cols_offsets[None, :]).ravel()
            data = np.tile([1.0, -2.0, 1.0], idx.size)
            fda = fda + sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()

        rhs = -fda[:, known_list].dot(a_flat[known_list])
        k = np.where(fda[:, nan_list].getnnz(axis=1) > 0)[0]
        A_nan = fda[k][:, nan_list]
        sol = lsqr(A_nan, rhs[k])[0]
        b_flat = a_flat.copy()
        b_flat[nan_list] = sol
        return b_flat.reshape((n, m), order="F")

    if method == 4:
        hv_list = np.array([[-1, -1, 0], [1, 1, 0], [-n, 0, -1], [n, 0, 1]])
        hv_springs = []
        for hv in hv_list:
            hvs = nan_list_rc + hv
            k = (hvs[:, 1] >= 0) & (hvs[:, 1] < n) & (hvs[:, 2] >= 0) & (hvs[:, 2] < m)
            hv_springs.append(np.column_stack([nan_list_rc[k, 0], hvs[k, 0]]))
        if hv_springs:
            hv_springs = np.vstack(hv_springs)
        else:
            hv_springs = np.zeros((0, 2), dtype=int)

        hv_springs = np.unique(np.sort(hv_springs, axis=1), axis=0)
        nhv = hv_springs.shape[0]
        rows_idx = np.repeat(np.arange(nhv), 2)
        cols_idx = hv_springs.ravel()
        data = np.tile([1.0, -1.0], nhv)
        springs = sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nhv, nm)).tocsr()

        rhs = -springs[:, known_list].dot(a_flat[known_list])
        A_nan = springs[:, nan_list]
        sol = lsqr(A_nan, rhs)[0]
        b_flat = a_flat.copy()
        b_flat[nan_list] = sol
        return b_flat.reshape((n, m), order="F")

    # method == 5
    fda = sparse.csr_matrix((nm, nm))
    # -1,-1
    L = (nan_list_rc[:, 1] > 0) & (nan_list_rc[:, 2] > 0)
    if np.any(L):
        idx = nan_list_rc[L, 0]
        rows_idx = np.repeat(idx, 2)
        cols_idx = np.vstack([idx - (n + 1), idx]).T.ravel()
        data = np.tile([1.0, -1.0], idx.size)
        fda = fda + sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()
    # 0,-1
    L = nan_list_rc[:, 2] > 0
    if np.any(L):
        idx = nan_list_rc[L, 0]
        rows_idx = np.repeat(idx, 2)
        cols_idx = np.vstack([idx - n, idx]).T.ravel()
        data = np.tile([1.0, -1.0], idx.size)
        fda = fda + sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()
    # +1,-1
    L = (nan_list_rc[:, 1] < (n - 1)) & (nan_list_rc[:, 2] > 0)
    if np.any(L):
        idx = nan_list_rc[L, 0]
        rows_idx = np.repeat(idx, 2)
        cols_idx = np.vstack([idx - (n - 1), idx]).T.ravel()
        data = np.tile([1.0, -1.0], idx.size)
        fda = fda + sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()
    # -1,0
    L = nan_list_rc[:, 1] > 0
    if np.any(L):
        idx = nan_list_rc[L, 0]
        rows_idx = np.repeat(idx, 2)
        cols_idx = np.vstack([idx - 1, idx]).T.ravel()
        data = np.tile([1.0, -1.0], idx.size)
        fda = fda + sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()
    # +1,0
    L = nan_list_rc[:, 1] < (n - 1)
    if np.any(L):
        idx = nan_list_rc[L, 0]
        rows_idx = np.repeat(idx, 2)
        cols_idx = np.vstack([idx + 1, idx]).T.ravel()
        data = np.tile([1.0, -1.0], idx.size)
        fda = fda + sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()
    # -1,+1
    L = (nan_list_rc[:, 1] > 0) & (nan_list_rc[:, 2] < (m - 1))
    if np.any(L):
        idx = nan_list_rc[L, 0]
        rows_idx = np.repeat(idx, 2)
        cols_idx = np.vstack([idx + (n - 1), idx]).T.ravel()
        data = np.tile([1.0, -1.0], idx.size)
        fda = fda + sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()
    # 0,+1
    L = nan_list_rc[:, 2] < (m - 1)
    if np.any(L):
        idx = nan_list_rc[L, 0]
        rows_idx = np.repeat(idx, 2)
        cols_idx = np.vstack([idx + n, idx]).T.ravel()
        data = np.tile([1.0, -1.0], idx.size)
        fda = fda + sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()
    # +1,+1
    L = (nan_list_rc[:, 1] < (n - 1)) & (nan_list_rc[:, 2] < (m - 1))
    if np.any(L):
        idx = nan_list_rc[L, 0]
        rows_idx = np.repeat(idx, 2)
        cols_idx = np.vstack([idx + (n + 1), idx]).T.ravel()
        data = np.tile([1.0, -1.0], idx.size)
        fda = fda + sparse.coo_matrix((data, (rows_idx, cols_idx)), shape=(nm, nm)).tocsr()

    rhs = -fda[:, known_list].dot(a_flat[known_list])
    b_flat = a_flat.copy()
    b_flat[nan_list] = spsolve(fda[nan_list][:, nan_list], rhs[nan_list])
    return b_flat.reshape((n, m), order="F")


def _identify_neighbors(n: int, m: int, nan_list_rc: np.ndarray, talks_to: np.ndarray) -> np.ndarray:
    """Identify neighbors of NaNs (0-based)."""
    if nan_list_rc.size == 0:
        return np.zeros((0, 3), dtype=int)

    nan_count = nan_list_rc.shape[0]
    talk_count = talks_to.shape[0]
    nn = np.zeros((nan_count * talk_count, 2), dtype=int)
    j0 = 0
    for i in range(talk_count):
        j1 = j0 + nan_count
        nn[j0:j1, :] = nan_list_rc[:, 1:3] + talks_to[i]
        j0 = j1

    L = (nn[:, 0] < 0) | (nn[:, 0] >= n) | (nn[:, 1] < 0) | (nn[:, 1] >= m)
    nn = nn[~L]
    lin = nn[:, 0] + nn[:, 1] * n
    neighbors_list = np.column_stack([lin, nn])
    neighbors_list = np.unique(neighbors_list, axis=0)
    # drop those which are also NaNs
    nan_set = set(map(tuple, nan_list_rc))
    neighbors_list = np.array([row for row in neighbors_list if tuple(row) not in nan_set], dtype=int)
    return neighbors_list
