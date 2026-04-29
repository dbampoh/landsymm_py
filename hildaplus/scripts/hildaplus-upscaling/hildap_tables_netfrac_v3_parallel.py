import argparse
import os
import sys
from multiprocessing import Pool
from pathlib import Path
from timeit import default_timer as timer

import numpy as np
from netCDF4 import Dataset

# Configuration - paths computed relative to script location (landsymm_py root is 3 levels up)
_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_HILDA_DATA = str(_PROJECT_ROOT / 'data' / 'geodata_py' / 'HILDA+' / 'data')
FNAME_DATAFILE = os.path.join(_HILDA_DATA, 'hildaplus_smoothed.nc')
FNAME_FM_DATAFILE = None  # Path to forest management file if separate (optional)
OUTPUT_DIR = os.path.join(_HILDA_DATA, 'output')
FNAME_GRIDLIST = os.path.join(_HILDA_DATA, 'input_data', 'gridlist_in_62892_and_climate.txt')

# ---------------------------------------------------------------------------
# YAML-driven HILDA+ → LPJ-GUESS mapping config
# ---------------------------------------------------------------------------
# See `hildaplus/config/README.md` for the schema. Default behavior is
# unchanged from the historical hardcoded mapping (parity verified by
# `landsymm/tests/test_landcover_config_parity.py`).
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
from hildaplus.config import (  # noqa: E402  (path-tweak-then-import is intentional)
    load_landcover_config,
    get_categories,
    get_forest_types,
    get_forest_management_codes,
    get_output_settings,
)

# Module-level placeholders populated by `_init_config()` once CLI args
# are resolved. These mirror the original hardcoded constants so existing
# in-script references keep working.
CFG = None
DLON_LPJG = 0.5
DLAT_LPJG = 0.5
PRECISION = 7
LC_FOREST_UNMANAGED = []
LC_FOREST_MANAGED = []
LC_GRASSLAND_SHRUBLAND = []
LC_SPARSE_BARREN = []
FOREST_TYPES = {}
FOREST_TYPES_MANAGED = {}
FOREST_TYPES_UNMANAGED = {}
# FM sentinel constants (previously undefined — latent bug fix).
FM_NO_FOREST = 0
FM_MANAGED = 1
FM_UNMANAGED = 2


def _init_config(config_path=None, profile=None):
    """Load the YAML mapping config and populate module-level shorthands."""
    global CFG, DLON_LPJG, DLAT_LPJG, PRECISION
    global LC_FOREST_UNMANAGED, LC_FOREST_MANAGED
    global LC_GRASSLAND_SHRUBLAND, LC_SPARSE_BARREN
    global FOREST_TYPES, FOREST_TYPES_MANAGED, FOREST_TYPES_UNMANAGED
    global FM_NO_FOREST, FM_MANAGED, FM_UNMANAGED

    CFG = load_landcover_config(config_path=config_path, profile=profile)

    out = get_output_settings(CFG)
    DLON_LPJG = out["dlon"]
    DLAT_LPJG = out["dlat"]
    PRECISION = out["precision"]

    raw = CFG["hilda_classes"]
    LC_FOREST_UNMANAGED = list(raw.get("forest_unmanaged", []))
    LC_FOREST_MANAGED = list(raw.get("forest_managed", []))
    LC_GRASSLAND_SHRUBLAND = list(raw.get("grassland_shrubland", []))
    LC_SPARSE_BARREN = list(raw.get("sparse_barren", []))

    FOREST_TYPES = get_forest_types(CFG, "combined")
    FOREST_TYPES_MANAGED = get_forest_types(CFG, "managed")
    FOREST_TYPES_UNMANAGED = get_forest_types(CFG, "unmanaged")

    fm_codes = get_forest_management_codes(CFG)
    FM_NO_FOREST = fm_codes.get("no_forest", 0)
    FM_MANAGED = fm_codes.get("managed", 1)
    FM_UNMANAGED = fm_codes.get("unmanaged", 2)


def get_lpjguess_categories(has_fm, forest_mode):
    """Backwards-compatible shim around `hildaplus.config.get_categories`."""
    if CFG is None:
        _init_config()
    return get_categories(CFG, has_fm=has_fm, forest_mode=forest_mode)


def read_gridlist(fname, start=0, nlines=-1):
    with open(fname) as f:
        gridlist = []
        for iline, line in enumerate(f):
            if iline < start:
                continue
            if iline == start + nlines:
                break
            bits = line.split()
            lon, lat = float(bits[0]), float(bits[1])
            gridlist.append({'lon': lon, 'lat': lat})
    return gridlist


def count_gridlist_lines(fname):
    count = 0
    with open(fname) as f:
        for line in f:
            if line.strip():
                count += 1
    return count


def get_dataset(fname, varname='LULC_states'):
    data = Dataset(fname)
    lon = data['longitude'][:]
    lat = data['latitude'][:]
    time = data['time'][:]
    dlon = lon[1] - lon[0]
    dlat = lat[1] - lat[0]
    return data[varname], lon, lat, time, dlon, dlat, data


def get_forest_management_data(fname):
    if fname and os.path.exists(fname):
        data = Dataset(fname)
        return data['forest_management'], data
    return None, None


def get_site_data(data_handle, lon, lat, dlon, dlat, site):
    lon_slice = np.where((lon >= site['lon'] - DLON_LPJG/2. + dlon/2.)
                        & (lon <= site['lon'] + DLON_LPJG/2. - dlon/2.))[0]
    lat_slice = np.where((lat >= site['lat'] - DLAT_LPJG/2. + dlat/2.)
                        & (lat <= site['lat'] + DLAT_LPJG/2. - dlat/2.))[0]
    data_site = (data_handle[:, lat_slice, lon_slice].data).astype('uint32')
    return data_site


def process_with_forest_management(lc_current, fm_current, categories):
    npixels = lc_current.size
    row = []
    for cat_name, lc_codes in categories.items():
        if cat_name == 'MANAGED_FOREST':
            forest_mask = np.isin(lc_current, LC_FOREST_UNMANAGED)
            managed_count = np.count_nonzero(forest_mask & (fm_current == FM_MANAGED))
            row.append(managed_count / npixels)
        elif cat_name == 'UNMANAGED_FOREST':
            forest_mask = np.isin(lc_current, LC_FOREST_UNMANAGED)
            unmanaged_count = np.count_nonzero(forest_mask & (fm_current == FM_UNMANAGED))
            row.append(unmanaged_count / npixels)
        elif cat_name == 'NATURAL':
            grass_sparse_count = np.count_nonzero(
                np.isin(lc_current, LC_GRASSLAND_SHRUBLAND + LC_SPARSE_BARREN)
            )
            forest_mask = np.isin(lc_current, LC_FOREST_UNMANAGED)
            unmanaged_forest_count = np.count_nonzero(forest_mask & (fm_current == FM_UNMANAGED))
            total_natural = grass_sparse_count + unmanaged_forest_count
            row.append(total_natural / npixels)
        else:
            n_occurrences = np.count_nonzero(np.isin(lc_current, lc_codes))
            row.append(n_occurrences / npixels)
    return row


def make_netfrac_table(gridlist, lc_data, lon, lat, time, dlon, dlat, fm_data=None,
                       integrated_fm=False, forest_mode="combined"):
    netfrac_table = []
    forestfrac_table = []
    categories = get_lpjguess_categories(has_fm=fm_data is not None, forest_mode=forest_mode)

    for site in gridlist:
        lc_data_site = get_site_data(lc_data, lon, lat, dlon, dlat, site)
        if fm_data is not None:
            fm_data_site = get_site_data(fm_data, lon, lat, dlon, dlat, site)

        npixels = lc_data_site[0].size
        for year_idx, year in enumerate(time[1:], 1):
            lc_current = lc_data_site[year_idx]
            row_netfrac = [site['lon'], site['lat'], year]

            for _, lc_codes in categories.items():
                n_occurrences = np.count_nonzero(np.isin(lc_current, lc_codes))
                row_netfrac.append(n_occurrences / npixels)

            netfrac_table.append(row_netfrac)

            forest_codes = LC_FOREST_UNMANAGED + (LC_FOREST_MANAGED if integrated_fm else [])
            forest_total = np.count_nonzero(np.isin(lc_current, forest_codes))
            row_forestfrac = [site['lon'], site['lat'], year]
            if forest_mode == "split" and integrated_fm:
                if forest_total == 0:
                    row_forestfrac.extend([0.0] * (len(FOREST_TYPES_MANAGED) + len(FOREST_TYPES_UNMANAGED)))
                else:
                    for _, codes in FOREST_TYPES_MANAGED.items():
                        count = np.count_nonzero(np.isin(lc_current, codes))
                        row_forestfrac.append(count / forest_total)
                    for _, codes in FOREST_TYPES_UNMANAGED.items():
                        count = np.count_nonzero(np.isin(lc_current, codes))
                        row_forestfrac.append(count / forest_total)
            else:
                if forest_total == 0:
                    row_forestfrac.extend([0.0] * len(FOREST_TYPES))
                else:
                    for _, codes in FOREST_TYPES.items():
                        count = np.count_nonzero(np.isin(lc_current, codes))
                        row_forestfrac.append(count / forest_total)
            forestfrac_table.append(row_forestfrac)

    return netfrac_table, forestfrac_table


def write_netfrac_table(table, fname, categories, include_header=True):
    width_coord = 10
    width_year = 6
    width_data = max(PRECISION + 5, max(len(c) for c in categories.keys()) + 1)

    cat_names = list(categories.keys())
    ndata_columns = len(cat_names)

    widths = [width_coord, width_coord, width_year] + [width_data] * ndata_columns
    format_list = [f'%{width_coord}.3f', f'%{width_coord}.3f', f'%{width_year}d'] \
                + [f'%{width_data}.{PRECISION}f'] * ndata_columns
    columns = ['lon', 'lat', 'year'] + cat_names

    header = ''
    for col, width in zip(columns, widths):
        header += col.rjust(width)

    table = np.array(table)
    table[:, 3:] = table[:, 3:].round(PRECISION)

    np.savetxt(fname, table, fmt=format_list, delimiter='', header=header if include_header else '',
               comments='' if include_header else '')


def process_chunk(args):
    (chunk_idx, gridlist, datafile, fmfile, integrated_fm, temp_dir,
     forest_mode, mapping_config, mapping_profile) = args

    # Workers may be spawned (not forked) on macOS/Windows, in which case
    # module-level config state isn't inherited from the parent. Re-init
    # safely either way.
    _init_config(config_path=mapping_config, profile=mapping_profile)

    lc_data, lon, lat, time, dlon, dlat, ds_lc = get_dataset(datafile, 'LULC_states')
    fm_data, ds_fm = get_forest_management_data(fmfile)

    netfrac_table, forestfrac_table = make_netfrac_table(
        gridlist, lc_data, lon, lat, time, dlon, dlat, fm_data, integrated_fm, forest_mode
    )

    if ds_fm:
        ds_fm.close()
    ds_lc.close()

    temp_path = os.path.join(temp_dir, f'netfrac_chunk_{chunk_idx:04d}.txt')
    categories = get_lpjguess_categories(has_fm=fm_data is not None, forest_mode=forest_mode)
    write_netfrac_table(netfrac_table, temp_path, categories, include_header=False)
    forest_path = os.path.join(temp_dir, f'forestfrac_chunk_{chunk_idx:04d}.txt')
    forest_headers = (
        {**FOREST_TYPES_MANAGED, **FOREST_TYPES_UNMANAGED}
        if forest_mode == "split" and integrated_fm
        else FOREST_TYPES
    )
    write_netfrac_table(forestfrac_table, forest_path, forest_headers, include_header=False)
    return chunk_idx, temp_path, forest_path


def parse_args():
    parser = argparse.ArgumentParser(description="Parallel netfrac upscaling (v3)")
    parser.add_argument("start", nargs="?", type=int, help="Gridlist start line (0-based)")
    parser.add_argument("nlines", nargs="?", type=int, help="Number of gridlist lines to process")
    parser.add_argument("suffix", nargs="?", default="", help="Output suffix")
    parser.add_argument("output_path", nargs="?", default=None, help="Override output file path")
    parser.add_argument("--benchmark-lines", type=int, default=0,
                        help="Process N gridlist lines, report timing, and exit (0 disables)")
    parser.add_argument("--workers", type=int, default=8, help="Number of worker processes")
    parser.add_argument("--chunk-lines", type=int, default=500,
                        help="Gridlist lines per worker chunk")
    parser.add_argument("--temp-dir", default=None,
                        help="Temp directory for chunk outputs (default: <output_dir>/tmp_netfrac_chunks)")
    parser.add_argument("--datafile", default=None, help="Path to smoothed LULC NetCDF")
    parser.add_argument("--gridlist", default=None, help="Path to gridlist file")
    parser.add_argument("--fmfile", default=None, help="Path to forest management NetCDF (optional)")
    parser.add_argument("--output", default=None, help="Override netfrac output file path (flag version)")
    parser.add_argument("--forestfrac-output", default=None,
                        help="Override forestfrac output file path (flag version)")
    parser.add_argument("--forest-mode", choices=["combined", "split"], default="combined",
                        help="Forest handling in outputs (default: combined)")
    parser.add_argument("--mapping-config", default=None,
                        help="Path to a custom HILDA+ -> LPJ-GUESS mapping YAML "
                             "(see hildaplus/config/README.md for the schema)")
    parser.add_argument("--mapping-profile", default=None,
                        help="Named mapping profile (e.g. lpjg_v3_default, lpjg_legacy_v1, "
                             "lpjg_treecrops_as_forest); resolves to "
                             "hildaplus/config/profiles/<name>.yaml")
    return parser.parse_args()


def main():
    args = parse_args()
    if (args.start is None) ^ (args.nlines is None):
        exit('Usage: python hildap_tables_netfrac_v3_parallel.py [start nlines [suffix [output_path]]] '
             '[--benchmark-lines N] [--workers N] [--chunk-lines N]')

    # Initialize the YAML-driven mapping config (CLI flags override env vars
    # which override the on-disk default; see hildaplus/config/README.md).
    _init_config(config_path=args.mapping_config, profile=args.mapping_profile)
    print(f"Land-cover mapping config: {CFG.get('_source_path')}")
    print(f"  Profile: {CFG.get('profile', {}).get('name', '<unnamed>')}")

    suffix = f'.{args.suffix}' if args.suffix != '' else ''
    output_path = args.output if args.output is not None else args.output_path
    if output_path:
        fname_netfrac = output_path
    else:
        fname_netfrac = os.path.join(OUTPUT_DIR, 'hildaplus_netfrac_1901_2020.txt' + suffix)

    temp_dir = args.temp_dir or os.path.join(os.path.dirname(fname_netfrac), 'tmp_netfrac_chunks')
    os.makedirs(temp_dir, exist_ok=True)

    datafile = args.datafile or FNAME_DATAFILE
    gridlist_path = args.gridlist or FNAME_GRIDLIST
    fmfile = args.fmfile if args.fmfile is not None else FNAME_FM_DATAFILE

    if args.start is None or args.nlines is None:
        gridlist = read_gridlist(gridlist_path)
    else:
        gridlist = read_gridlist(gridlist_path, start=args.start, nlines=args.nlines)

    total_gridlist_lines = None
    if args.benchmark_lines and args.benchmark_lines > 0:
        total_gridlist_lines = count_gridlist_lines(gridlist_path)
        gridlist = gridlist[:args.benchmark_lines]

    # Determine integrated FM
    lc_data, lon, lat, time, dlon, dlat, ds_lc = get_dataset(datafile, 'LULC_states')
    sample_data = get_site_data(lc_data, lon, lat, dlon, dlat, gridlist[0])
    integrated_fm = np.any(np.isin(sample_data, LC_FOREST_MANAGED))
    ds_lc.close()
    forest_mode = args.forest_mode
    if forest_mode == "split" and not integrated_fm:
        print("Warning: forest_mode=split requested but no integrated FM codes found. Using combined forest.")
        forest_mode = "combined"

    # Build chunks
    # NOTE: We pass the resolved config path through so worker processes can
    # re-load the YAML deterministically. (Sending the parsed dict via pickle
    # would also work but introduces unnecessary cross-process serialization.)
    chunks = []
    resolved_cfg_path = CFG.get('_source_path') if CFG else None
    for idx in range(0, len(gridlist), args.chunk_lines):
        chunk_gridlist = gridlist[idx:idx + args.chunk_lines]
        chunk_idx = idx // args.chunk_lines
        chunks.append((chunk_idx, chunk_gridlist, datafile, fmfile, integrated_fm,
                       temp_dir, forest_mode, resolved_cfg_path, None))

    categories = get_lpjguess_categories(has_fm=fmfile is not None, forest_mode=forest_mode)
    start_time = timer()

    with Pool(processes=args.workers) as pool:
        results = list(pool.imap_unordered(process_chunk, chunks))

    # Merge chunk outputs in order
    results.sort(key=lambda x: x[0])
    with open(fname_netfrac, 'w') as out_f:
        # Write header only (no dummy data row)
        width_coord = 10
        width_year = 6
        width_data = max(PRECISION + 5, max(len(c) for c in categories.keys()) + 1)
        columns = ['lon', 'lat', 'year'] + list(categories.keys())
        widths = [width_coord, width_coord, width_year] + [width_data] * len(categories)
        header = ''.join([col.rjust(width) for col, width in zip(columns, widths)])
        out_f.write(header + "\n")

        # Append chunk files
        for _, temp_path, _ in results:
            with open(temp_path, 'r') as tf:
                out_f.write(tf.read())
            os.remove(temp_path)

    if args.forestfrac_output:
        fname_forestfrac = args.forestfrac_output
    else:
        fname_forestfrac = os.path.join(os.path.dirname(fname_netfrac), 'hildaplus_forestfrac_1901_2020.txt' + suffix)
    with open(fname_forestfrac, 'w') as out_f:
        width_coord = 10
        width_year = 6
        forest_headers = (
            {**FOREST_TYPES_MANAGED, **FOREST_TYPES_UNMANAGED}
            if forest_mode == "split" and integrated_fm
            else FOREST_TYPES
        )
        width_data = max(PRECISION + 5, max(len(c) for c in forest_headers.keys()) + 1)
        columns = ['lon', 'lat', 'year'] + list(forest_headers.keys())
        widths = [width_coord, width_coord, width_year] + [width_data] * len(forest_headers)
        header = ''.join([col.rjust(width) for col, width in zip(columns, widths)])
        out_f.write(header + "\n")

        for _, _, forest_path in results:
            with open(forest_path, 'r') as tf:
                out_f.write(tf.read())
            os.remove(forest_path)

    elapsed = timer() - start_time

    if args.benchmark_lines and args.benchmark_lines > 0 and total_gridlist_lines:
        per_line = elapsed / len(gridlist)
        est_total = per_line * total_gridlist_lines
        print(f"Benchmark complete: {len(gridlist)} lines in {elapsed:.1f}s "
              f"({per_line:.4f}s/line). Estimated full run: {est_total/3600:.2f} hours.")
        return

    print(f"Wrote output to {fname_netfrac} in {elapsed:.1f}s")
    print(f"Wrote output to {fname_forestfrac} in {elapsed:.1f}s")


if __name__ == '__main__':
    main()
