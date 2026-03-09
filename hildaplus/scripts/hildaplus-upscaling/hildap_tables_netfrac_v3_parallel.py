import argparse
import os
from multiprocessing import Pool
from pathlib import Path
from timeit import default_timer as timer

import numpy as np
from netCDF4 import Dataset

# Configuration - paths computed relative to script location (landsymm_py root is 3 levels up)
_HILDA_DATA = str(Path(__file__).resolve().parents[3] / 'data' / 'geodata_py' / 'HILDA+' / 'data')
FNAME_DATAFILE = os.path.join(_HILDA_DATA, 'hildaplus_smoothed.nc')
FNAME_FM_DATAFILE = None  # Path to forest management file if separate (optional)
DLON_LPJG = 0.5  # LPJ-GUESS resolution in longitude (deg)
DLAT_LPJG = 0.5  # LPJ-GUESS resolution in latitude (deg)
PRECISION = 7    # Number of decimals to round output
OUTPUT_DIR = os.path.join(_HILDA_DATA, 'output')
FNAME_GRIDLIST = os.path.join(_HILDA_DATA, 'input_data', 'gridlist_in_62892_and_climate.txt')

# Land cover type definitions for new HILDA+ v2
LC_OCEAN = [0]
LC_URBAN = [11]
LC_ANNUAL_CROPS = [22]
LC_TREE_CROPS = [23]
LC_AGROFORESTRY = [24]
LC_PASTURE = [33]
LC_FOREST_UNMANAGED = [40, 41, 42, 43, 44, 45]  # Natural forest codes
LC_FOREST_MANAGED = [400, 410, 420, 430, 440, 450]  # Managed forest codes (×10)
LC_GRASSLAND_SHRUBLAND = [55]
LC_SPARSE_BARREN = [66]
LC_WATER = [77]
LC_NODATA = [99]

# Forest fractions by type (codes include integrated managed variants if present)
FOREST_TYPES = {
    "ForestNE": [41, 410],
    "ForestND": [43, 430],
    "ForestBE": [42, 420],
    "ForestBD": [44, 440],
    "ForestPNV": [40, 45, 400, 450],
}

FOREST_TYPES_MANAGED = {
    "ForestNE_MANAGED": [410],
    "ForestND_MANAGED": [430],
    "ForestBE_MANAGED": [420],
    "ForestBD_MANAGED": [440],
    "ForestPNV_MANAGED": [400, 450],
}

FOREST_TYPES_UNMANAGED = {
    "ForestNE_UNMANAGED": [41],
    "ForestND_UNMANAGED": [43],
    "ForestBE_UNMANAGED": [42],
    "ForestBD_UNMANAGED": [44],
    "ForestPNV_UNMANAGED": [40, 45],
}


def get_lpjguess_categories(has_fm, forest_mode):
    categories = {
        'URBAN': LC_URBAN,
        'CROPLAND': LC_ANNUAL_CROPS + LC_TREE_CROPS + LC_AGROFORESTRY,
        'PASTURE': LC_PASTURE,
    }
    if forest_mode == "split" and has_fm:
        categories.update({
            'FOREST_MANAGED': LC_FOREST_MANAGED,
            'FOREST_UNMANAGED': LC_FOREST_UNMANAGED,
            'NATURAL': LC_GRASSLAND_SHRUBLAND + LC_SPARSE_BARREN,
        })
    else:
        categories.update({
            'FOREST': LC_FOREST_UNMANAGED + LC_FOREST_MANAGED,
            'NATURAL': LC_GRASSLAND_SHRUBLAND + LC_SPARSE_BARREN,
        })
    categories.update({
        'BARREN': LC_OCEAN + LC_WATER + LC_NODATA,
    })
    return categories


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
    chunk_idx, gridlist, datafile, fmfile, integrated_fm, temp_dir, forest_mode = args

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
    return parser.parse_args()


def main():
    args = parse_args()
    if (args.start is None) ^ (args.nlines is None):
        exit('Usage: python hildap_tables_netfrac_v3_parallel.py [start nlines [suffix [output_path]]] '
             '[--benchmark-lines N] [--workers N] [--chunk-lines N]')

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
    chunks = []
    for idx in range(0, len(gridlist), args.chunk_lines):
        chunk_gridlist = gridlist[idx:idx + args.chunk_lines]
        chunk_idx = idx // args.chunk_lines
        chunks.append((chunk_idx, chunk_gridlist, datafile, fmfile, integrated_fm, temp_dir, forest_mode))

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
