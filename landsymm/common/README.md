# `landsymm.common` — Shared Utilities

Shared utility modules used across all pipeline stages.

## Modules

| Module | Description |
|--------|-------------|
| `aggregation.py` | Spatial aggregation (e.g., 5-arcmin → 0.5° grid) |
| `geoarray.py` | `GeoArray` dataclass for georeferenced arrays |
| `gridlist.py` | Load and filter LPJ-GUESS gridlist files |
| `interpolation.py` | Inpainting / gap-filling for gridded data |
| `io_utils.py` | File I/O helpers (delete existing outputs, etc.) |
| `logging.py` | Standardized logging setup |
| `lpjg_io.py` | Core I/O for LPJ-GUESS table/map formats and `.mat` files |
| `mapping_tools.py` | Coordinate transformations (yxz ↔ xz, vector ↔ map) |
| `validation.py` | Data validation utilities |

## Usage

```python
from landsymm.common.lpjg_io import read_table, save_table, maps2table
from landsymm.common.aggregation import aggregate
from landsymm.common.interpolation import inpaint_nans
from landsymm.common.geoarray import GeoArray
from landsymm.common.mapping_tools import yxz_to_xz, xz_to_yxz
```
