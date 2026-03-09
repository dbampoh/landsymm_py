import argparse
import numpy as np
from netCDF4 import Dataset
from pathlib import Path
from timeit import default_timer as timer
from sys import argv, version_info
import os

# Configuration - paths computed relative to script location (landsymm_py root is 3 levels up)
# This script is a netfrac-only adaptation of the original upscaling workflow.
_HILDA_DATA = str(Path(__file__).resolve().parents[3] / 'data' / 'geodata_py' / 'HILDA+' / 'data')
FNAME_DATAFILE = os.path.join(_HILDA_DATA, 'hildaplus_smoothed.nc')
FNAME_FM_DATAFILE = None  # Path to forest management file if separate (optional)
DLON_LPJG = 0.5  # LPJ-GUESS resolution in longitude (deg)
DLAT_LPJG = 0.5  # LPJ-GUESS resolution in latitude (deg)
PRECISION = 7    # Number of decimals to round output
OUTPUT_DIR = os.path.join(_HILDA_DATA, 'output')

# Gridlist file
FNAME_GRIDLIST = os.path.join(_HILDA_DATA, 'input_data', 'gridlist_in_62892_and_climate.txt')

# Land cover type definitions for new HILDA+ v2
# Individual land cover codes
LC_OCEAN = [0]
LC_URBAN = [11]
LC_ANNUAL_CROPS = [22]
LC_TREE_CROPS = [23]
LC_AGROFORESTRY = [24]
LC_PASTURE = [33]
LC_FOREST_UNMANAGED = [40, 41, 42, 43, 44, 45]  # Forest codes
LC_FOREST_MANAGED = [400, 410, 420, 430, 440, 450]  # Managed forest codes (×10, if integrated)
LC_GRASSLAND_SHRUBLAND = [55]
LC_SPARSE_BARREN = [66]
LC_WATER = [77]
LC_NODATA = [99]

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

# Aggregated categories for LPJ-GUESS (final mapping agreed in this workflow).
def get_lpjguess_categories(has_fm, forest_mode):
    """Define LPJ-GUESS land cover categories based on new HILDA+ codes."""
    categories = {
        'URBAN': LC_URBAN,
        'CROPLAND': LC_ANNUAL_CROPS + LC_TREE_CROPS + LC_AGROFORESTRY,
        'PASTURE': LC_PASTURE,
    }
    if forest_mode == "split" and has_fm:
        categories.update({
            "FOREST_MANAGED": LC_FOREST_MANAGED,
            "FOREST_UNMANAGED": LC_FOREST_UNMANAGED,
            "NATURAL": LC_GRASSLAND_SHRUBLAND + LC_SPARSE_BARREN,
        })
    else:
        categories.update({
            "FOREST": LC_FOREST_UNMANAGED + LC_FOREST_MANAGED,
            "NATURAL": LC_GRASSLAND_SHRUBLAND + LC_SPARSE_BARREN,
        })
    categories.update({
        'BARREN': LC_OCEAN + LC_WATER + LC_NODATA,  # Only truly barren areas
    })
    return categories


def read_gridlist(fname, start=0, nlines=-1):
    """Read gridlist file containing lon/lat coordinates."""
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
    """Count non-empty lines in a gridlist file."""
    count = 0
    with open(fname) as f:
        for line in f:
            if line.strip():
                count += 1
    return count


def get_dataset(fname, varname='LULC_states'):
    """Open NetCDF dataset and extract coordinates."""
    data = Dataset(fname)
    lon = data['longitude'][:]
    lat = data['latitude'][:]
    time = data['time'][:]
    dlon = lon[1] - lon[0]
    dlat = lat[1] - lat[0]
    return data[varname], lon, lat, time, dlon, dlat


def get_forest_management_data(fname):
    """Load forest management data if available."""
    if fname and os.path.exists(fname):
        data = Dataset(fname)
        return data['forest_management']
    return None


def get_site_data(data_handle, lon, lat, dlon, dlat, site):
    """Extract data for a specific site from the dataset."""
    lon_slice = np.where((lon >= site['lon'] - DLON_LPJG/2. + dlon/2.)
                        & (lon <= site['lon'] + DLON_LPJG/2. - dlon/2.))[0]
    lat_slice = np.where((lat >= site['lat'] - DLAT_LPJG/2. + dlat/2.)
                        & (lat <= site['lat'] + DLAT_LPJG/2. - dlat/2.))[0]
    data_site = (data_handle[:, lat_slice, lon_slice].data).astype('uint32')
    return data_site


def make_netfrac_table(gridlist, lc_data, lon, lat, time, dlon, dlat, fm_data=None, forest_mode="combined"):
    """Generate net fractions and forest fraction tables."""
    netfrac_table = []
    forestfrac_table = []
    categories = get_lpjguess_categories(has_fm=fm_data is not None, forest_mode=forest_mode)
    
    # Check if we have integrated forest management (codes 400-450 in LULC data)
    integrated_fm = False
    if fm_data is None:
        print("Checking for integrated forest management codes...")
        sample_data = get_site_data(lc_data, lon, lat, dlon, dlat, gridlist[0])
        if np.any(np.isin(sample_data, LC_FOREST_MANAGED)):
            integrated_fm = True
            print("Found integrated forest management (codes 400-450)")
        elif forest_mode == "split":
            print("Warning: forest_mode=split requested but no integrated FM codes found. Using combined forest.")
            forest_mode = "combined"
    
    for i, site in enumerate(gridlist):
        print(f"Processing site {i+1}/{len(gridlist)}: lon={site['lon']}, lat={site['lat']}")
        
        # Retrieve slice of data for this site
        lc_data_site = get_site_data(lc_data, lon, lat, dlon, dlat, site)
        
        # Get forest management data if available
        if fm_data is not None:
            fm_data_site = get_site_data(fm_data, lon, lat, dlon, dlat, site)
        
        npixels = lc_data_site[0].size
        
        # Process each year
        for year_idx, year in enumerate(time[1:], 1):
            lc_current = lc_data_site[year_idx]
            
            # Start row with coordinates and year
            row_netfrac = [site['lon'], site['lat'], year]
            
            # Netfrac (no FM-specific outputs)
            for _, lc_codes in categories.items():
                n_occurrences = np.count_nonzero(np.isin(lc_current, lc_codes))
                row_netfrac.append(n_occurrences / npixels)
            
            netfrac_table.append(row_netfrac)

            # Forest fractions (relative to forest only)
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


def write_netfrac_table(table, fname, categories):
    """Write the net fractions table to file."""
    width_coord = 10
    width_year = 6
    width_data = max(PRECISION + 5, max(len(c) for c in cat_names) + 1)
    
    cat_names = list(categories.keys())
    ndata_columns = len(cat_names)
    
    widths = [width_coord, width_coord, width_year] + [width_data] * ndata_columns
    format_list = [f'%{width_coord}.3f', f'%{width_coord}.3f', f'%{width_year}d'] \
                + [f'%{width_data}.{PRECISION}f'] * ndata_columns
    columns = ['lon', 'lat', 'year'] + cat_names
    
    # Create header
    header = ''
    for col, width in zip(columns, widths):
        header += col.rjust(width)
    
    # Convert to numpy array and round
    table = np.array(table)
    table[:, 3:] = table[:, 3:].round(PRECISION)
    
    # Save to file
    np.savetxt(fname, table, fmt=format_list, delimiter='', header=header, comments='')


def main(fname_dataset, fname_gridlist, fname_fm=None, start=None, nlines=None, suffix='',
         output_path=None, forestfrac_path=None, benchmark_lines=0, forest_mode="combined"):
    """Main function to generate netfrac file."""
    start_time = timer()
    # Output filename
    suffix = f'.{suffix}' if suffix != '' else ''
    if output_path:
        fname_netfrac = output_path
    else:
        fname_netfrac = os.path.join(OUTPUT_DIR, 'hildaplus_netfrac_1901_2020.txt' + suffix)
    if forestfrac_path:
        fname_forestfrac = forestfrac_path
    else:
        fname_forestfrac = os.path.join(OUTPUT_DIR, 'hildaplus_forestfrac_1901_2020.txt' + suffix)
    
    # Print runtime information
    if start is None:
        print(f"Reading all entries in gridlist.")
    else:
        print(f"Reading gridlist range [{start}:{start+nlines}].")
    print(f"Net fractions table will be written to {fname_netfrac}")
    print(f"Using new HILDA+ v2 land cover categories with separate crop types")
    if fname_fm:
        print(f"Forest management data: {fname_fm}")
    else:
        print("Forest management will be detected from integrated codes (400-450)")
    
    # Read LPJ-GUESS gridlist
    if start is None or nlines is None:
        gridlist = read_gridlist(fname_gridlist)
    else:
        gridlist = read_gridlist(fname_gridlist, start=start, nlines=nlines)

    total_gridlist_lines = None
    if benchmark_lines and benchmark_lines > 0:
        total_gridlist_lines = count_gridlist_lines(fname_gridlist)
        gridlist = gridlist[:benchmark_lines]
    
    print(f"Found {len(gridlist)} grid points to process")
    
    # Retrieve dataset
    print("Opening dataset...")
    data, lon, lat, time, dlon, dlat = get_dataset(fname_dataset, 'LULC_states')
    print(f"Dataset covers years: {time[0]} to {time[-1]}")
    
    # Load forest management data if available (separate FM file case).
    fm_data = None
    if fname_fm:
        print("Loading forest management data...")
        fm_data = get_forest_management_data(fname_fm)
    
    # Get categories
    categories = get_lpjguess_categories(has_fm=fm_data is not None, forest_mode=forest_mode)
    print("\nLand cover categories:")
    for cat, codes in categories.items():
        print(f"  {cat}: codes {codes}")
    
    # Construct netfrac table
    print("\nConstructing net fractions table...")
    netfrac_table, forestfrac_table = make_netfrac_table(
        gridlist, data, lon, lat, time, dlon, dlat, fm_data, forest_mode=forest_mode
    )
    
    # Write output
    print(f"Writing to {fname_netfrac}...")
    write_netfrac_table(netfrac_table, fname_netfrac, categories)
    print(f"Writing to {fname_forestfrac}...")
    if forest_mode == "split" and (fm_data is not None or integrated_fm):
        forest_headers = {**FOREST_TYPES_MANAGED, **FOREST_TYPES_UNMANAGED}
        write_netfrac_table(forestfrac_table, fname_forestfrac, forest_headers)
    else:
        write_netfrac_table(forestfrac_table, fname_forestfrac, FOREST_TYPES)

    if benchmark_lines and benchmark_lines > 0 and total_gridlist_lines:
        elapsed = timer() - start_time
        per_line = elapsed / len(gridlist) if len(gridlist) else 0
        est_total = per_line * total_gridlist_lines
        print(f"Benchmark complete: {len(gridlist)} lines in {elapsed:.1f}s "
              f"({per_line:.4f}s/line). Estimated full run: {est_total/3600:.2f} hours.")
        return

    print("Done!")


def parse_args():
    parser = argparse.ArgumentParser(description="HILDA+ netfrac upscaling (v3)")
    parser.add_argument("start", nargs="?", type=int, help="Gridlist start line (0-based)")
    parser.add_argument("nlines", nargs="?", type=int, help="Number of gridlist lines to process")
    parser.add_argument("suffix", nargs="?", default="", help="Output suffix")
    parser.add_argument("output_path", nargs="?", default=None, help="Override output file path")
    parser.add_argument("--benchmark-lines", type=int, default=0,
                        help="Process N gridlist lines, report timing, and exit (0 disables)")
    parser.add_argument("--datafile", default=None, help="Path to smoothed LULC NetCDF")
    parser.add_argument("--gridlist", default=None, help="Path to gridlist file")
    parser.add_argument("--fmfile", default=None, help="Path to forest management NetCDF (optional)")
    parser.add_argument("--output", default=None, help="Override netfrac output file path (flag version)")
    parser.add_argument("--forestfrac-output", default=None,
                        help="Override forestfrac output file path (flag version)")
    parser.add_argument("--forest-mode", choices=["combined", "split"], default="combined",
                        help="Forest handling in outputs (default: combined)")
    return parser.parse_args()


if __name__ == '__main__':
    # Check Python version
    if version_info < (3, 7):
        exit('Python version must be 3.7 or newer.')

    args = parse_args()
    if (args.start is None) ^ (args.nlines is None):
        exit('Usage: python hildap_tables_netfrac_v3.py [start nlines [suffix [output_path]]] [--benchmark-lines N]')

    # Start timer
    t_start = timer()

    # Run main function
    datafile = args.datafile or FNAME_DATAFILE
    gridlist = args.gridlist or FNAME_GRIDLIST
    fmfile = args.fmfile if args.fmfile is not None else FNAME_FM_DATAFILE
    output_path = args.output if args.output is not None else args.output_path
    forestfrac_path = args.forestfrac_output

    main(
        datafile,
        gridlist,
        fmfile,
        args.start,
        args.nlines,
        args.suffix,
        output_path,
        forestfrac_path,
        benchmark_lines=args.benchmark_lines,
        forest_mode=args.forest_mode,
    )

    # Print elapsed time
    t_end = timer()
    print(f'Elapsed time: {round(t_end - t_start, 3)} seconds')