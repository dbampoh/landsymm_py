"""Diagnostics and plotting (stub).

Pseudocode:
1) map_missing:
   - find NaNs across variables
   - convert to map (gridlist list2map)
   - render shaded map
2) shade_map:
   - pcolor + colorbar
   - save to file
"""

# TODO:
# - implement missing-cell maps and shaded diagnostics
# - ensure figure naming matches MATLAB outputs
from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np

from landsymm.common.mapping_tools import vector_to_map


def map_missing(garr, gridlist, title: str, fig_dir: str):
    """Generate missing-cell diagnostics (TODO)."""
    isbad_x = np.isnan(garr)
    while np.count_nonzero(np.array(isbad_x.shape) > 1) > 1:
        isbad_x = np.any(isbad_x, axis=1)
    isbad_yx = vector_to_map(isbad_x, gridlist["mask_YX"].shape, gridlist["list2map"])

    os.makedirs(fig_dir, exist_ok=True)
    fig_outfile = os.path.join(fig_dir, f"cells_missing_from_{title}.png")
    n_missing = int(np.sum(isbad_x))
    title = f"{n_missing} cells missing from {title}"
    shade_map(isbad_yx, title, fig_outfile)
    return isbad_yx


def shade_map(map_yx, title: str, outfile: str):
    """Plot shaded map for diagnostics (TODO)."""
    fig, ax = plt.subplots(figsize=(10, 6))
    cax = ax.imshow(map_yx, origin="lower", cmap="cool", vmin=0, vmax=1)
    ax.set_axis_off()
    ax.set_title(title.replace("_", "\\_"))
    fig.colorbar(cax, ax=ax, fraction=0.046, pad=0.04)
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)
