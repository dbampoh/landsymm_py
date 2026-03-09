import numpy as np
from netCDF4 import Dataset
from timeit import default_timer as timer
from sys import argv, version_info

FNAME_DATAFILE = '/home/belda-d/data/hilda_plus/smoothed/hildaplus_smoothed.nc'
DLON_LPJG = 0.5 # LPJ-GUESS resolution in longitude (deg)
DLAT_LPJG = 0.5 # LPJ-GUESS resolution in latitude (deg)

PRECISION = 7 # Number of decimals to round output
PRECISION = max(PRECISION, 6) # At least 6 decimals

# Gridlist file
FNAME_GRIDLIST = 'lu_1700_luh2_aggregate_sum2x2_midpoint_nourban_gcp2022_2022_08_10.txt.NOHEADER'
# Check that all net fractions add up to 1
NETFRAC_CHECK = True
# Check that all forest fractions add up to 1
FORESTFRAC_CHECK = True

# List-of-strings-to-array-uint32
def ls2a_u32(l):
    il = [int(x) for x in l]
    return np.array(il, dtype='uint32')


# Land cover types
# ----------------

U = ['11', ]                                # Urban
C = ['22', ]                                # Cropland
P = ['33', ]                                # Grass/shrubland
F = ['40', '41', '42', '43', '44', '45', ]  # Forest
G = ['55', ]                                # Grassland
B = ['66', ]                                # Barren


# Transition types
# ----------------

# | Transition | Code       |
# |------------|------------|
# | yx         | YX         |
# | xs         | XF, XG     |
# | xf         | XF0        |
# | sx         | FX, GX     |
# | fx         | F0X        |
# | vx         | F1X, G1X   |
# | sf         | FF0, GF0   |
# | fs         | F0F, F0G   |
# | vf         | F1F0, G1F0 |

# yx types, where y, x are any of U, C, P, B
# ------------------------------------------

UU = ls2a_u32([e1+e2 for e1 in U for e2 in U])
UC = ls2a_u32([e1+e2 for e1 in U for e2 in C])
UP = ls2a_u32([e1+e2 for e1 in U for e2 in P])
UB = ls2a_u32([e1+e2 for e1 in U for e2 in B])

CU = ls2a_u32([e1+e2 for e1 in C for e2 in U])
CC = ls2a_u32([e1+e2 for e1 in C for e2 in C])
CP = ls2a_u32([e1+e2 for e1 in C for e2 in P])
CB = ls2a_u32([e1+e2 for e1 in C for e2 in B])

PU = ls2a_u32([e1+e2 for e1 in P for e2 in U])
PC = ls2a_u32([e1+e2 for e1 in P for e2 in C])
PP = ls2a_u32([e1+e2 for e1 in P for e2 in P])
PB = ls2a_u32([e1+e2 for e1 in P for e2 in B])

BU = ls2a_u32([e1+e2 for e1 in B for e2 in U])
BC = ls2a_u32([e1+e2 for e1 in B for e2 in C])
BP = ls2a_u32([e1+e2 for e1 in B for e2 in P])
BB = ls2a_u32([e1+e2 for e1 in B for e2 in B])

# xs types
# --------

UF = ls2a_u32([e1+e2 for e1 in U for e2 in F])
CF = ls2a_u32([e1+e2 for e1 in C for e2 in F])
PF = ls2a_u32([e1+e2 for e1 in P for e2 in F])
BF = ls2a_u32([e1+e2 for e1 in B for e2 in F])

UG = ls2a_u32([e1+e2 for e1 in U for e2 in G])
CG = ls2a_u32([e1+e2 for e1 in C for e2 in G])
PG = ls2a_u32([e1+e2 for e1 in P for e2 in G])
BG = ls2a_u32([e1+e2 for e1 in B for e2 in G])

US = np.concatenate([UF, UG])
CS = np.concatenate([CF, CG])
PS = np.concatenate([PF, PG])
BS = np.concatenate([BF, BG])

# xf types
# --------

UF0 = ls2a_u32([e1+e2+'0' for e1 in U for e2 in F])
CF0 = ls2a_u32([e1+e2+'0' for e1 in C for e2 in F])
PF0 = ls2a_u32([e1+e2+'0' for e1 in P for e2 in F])
BF0 = ls2a_u32([e1+e2+'0' for e1 in B for e2 in F])

UF = UF0
CF = CF0
PF = PF0
BF = BF0

# sx types
# --------

FU = ls2a_u32([e1+e2 for e1 in F for e2 in U])
FC = ls2a_u32([e1+e2 for e1 in F for e2 in C])
FP = ls2a_u32([e1+e2 for e1 in F for e2 in P])
FB = ls2a_u32([e1+e2 for e1 in F for e2 in B])

GU  = ls2a_u32([e1+e2 for e1 in G for e2 in U])
GC  = ls2a_u32([e1+e2 for e1 in G for e2 in C])
GP  = ls2a_u32([e1+e2 for e1 in G for e2 in P])
GB  = ls2a_u32([e1+e2 for e1 in G for e2 in B])

SU = np.concatenate([FU, GU])
SC = np.concatenate([FC, GC])
SP = np.concatenate([FP, GP])
SB = np.concatenate([FB, GB])

# fx types
# --------

F0U = ls2a_u32([e1+'0'+e2 for e1 in F for e2 in U])
F0C = ls2a_u32([e1+'0'+e2 for e1 in F for e2 in C])
F0P = ls2a_u32([e1+'0'+e2 for e1 in F for e2 in P])
F0B = ls2a_u32([e1+'0'+e2 for e1 in F for e2 in B])

FU = F0U
FC = F0C
FP = F0P
FB = F0B

# vx types
# --------

F1U = ls2a_u32([e1+'1'+e2 for e1 in F for e2 in U])
F1C = ls2a_u32([e1+'1'+e2 for e1 in F for e2 in C])
F1P = ls2a_u32([e1+'1'+e2 for e1 in F for e2 in P])
F1B = ls2a_u32([e1+'1'+e2 for e1 in F for e2 in B])

G1U = ls2a_u32([e1+'1'+e2 for e1 in G for e2 in U])
G1C = ls2a_u32([e1+'1'+e2 for e1 in G for e2 in C])
G1P = ls2a_u32([e1+'1'+e2 for e1 in G for e2 in P])
G1B = ls2a_u32([e1+'1'+e2 for e1 in G for e2 in B])

VU = np.concatenate([F1U, G1U])
VC = np.concatenate([F1C, G1C])
VP = np.concatenate([F1P, G1P])
VB = np.concatenate([F1B, G1B])

# sf types
# --------

FF0 = ls2a_u32([e1+e2+'0' for e1 in F for e2 in F])
GF0 = ls2a_u32([e1+e2+'0' for e1 in G for e2 in F])

SF = np.concatenate([FF0, GF0])

# fs types
# --------

F0F = ls2a_u32([e1+'0'+e2 for e1 in F for e2 in F])
F0G = ls2a_u32([e1+'0'+e2 for e1 in F for e2 in G])

FS = np.concatenate([F0F, F0G])

# vf types
# --------

F1F0 = ls2a_u32([e1+'1'+e2+'0' for e1 in F for e2 in F])
G1F0 = ls2a_u32([e1+'1'+e2+'0' for e1 in G for e2 in F])

VF = np.concatenate([F1F0, G1F0])

# Gridcells that stay primary
# ---------------------------

F1F = ls2a_u32([e1+'1'+e2 for e1 in F for e2 in F])
F1G = ls2a_u32([e1+'1'+e2 for e1 in F for e2 in G])
G1F = ls2a_u32([e1+'1'+e2 for e1 in G for e2 in F])
G1G = ls2a_u32([e1+'1'+e2 for e1 in G for e2 in G])

STAYS_PRIMARY = np.concatenate([F1F, F1G, G1F, G1F, ])


# TRANSITION TYPE EQUIVALENCES
# ----------------------------

TRANSITIONS = [UP, UC, UB, PU, PC, PB, CU,
               CP, CB, BU, BP, BC, US, PS,
               CS, BS, UF, PF, CF, BF, SU,
               SP, SC, SB, FU, FP, FC, FB,
               VU, VP, VC, VB, SF, FS, VF,
               ]
TR_KEYS = ['up', 'uc', 'ub', 'pu', 'pc', 'pb', 'cu',
           'cp', 'cb', 'bu', 'bp', 'bc', 'us', 'ps',
           'cs', 'bs', 'uf', 'pf', 'cf', 'bf', 'su',
           'sp', 'sc', 'sb', 'fu', 'fp', 'fc', 'fb',
           'vu', 'vp', 'vc', 'vb', 'sf', 'fs', 'vf',
           ]
TRANSITIONS = dict(zip(TR_KEYS, TRANSITIONS))


# LAND COVER TYPE EQUIVALENCES
# ----------------------------

# | Code    | LPJ-GUESS LCT                           |                        
# |---------|-----------------------------------------|                        
# | 0       | Barren (ocean)                          |                        
# | 11      | Urban                                   |                        
# | 22      | Crop                                    |                        
# | 33      | Pasture                                 |                        
# | 40-45   | Natural (unmanaged forest)              |
# | 400-450 | Forest (managed forest)                 |
# | 55      | Natural                                 |                        
# | 66      | Barren (barren/other land)              |                        
# | 77      | Barren (water)                          |                        
# | 99      | Barren (no data)                        |

URBAN = ls2a_u32(U)
CROPLAND = ls2a_u32(C)
PASTURE = ls2a_u32(P)
FOREST = ls2a_u32(F)*10
NATURAL = np.concatenate([ls2a_u32(F), ls2a_u32(G)])
BARREN = ls2a_u32(B + ['77', '99', '0', ])

FOREST_TYPES = [[410, ], [430, ], [420, ], [440, ], [450, 400, ], ]
FOREST_KEYS = ['ForestNE', 'ForestND', 'ForestBE', 'ForestBD', 'ForestPNV', ]
FOREST_TYPES = dict(zip(FOREST_KEYS, FOREST_TYPES))
NFOREST_TYPES = len(FOREST_TYPES)

LC_TYPES = [URBAN, CROPLAND, PASTURE, FOREST, NATURAL, BARREN, ]
LC_KEYS = ['URBAN', 'CROPLAND', 'PASTURE', 'FOREST', 'NATURAL', 'BARREN', ]
LC_TYPES = dict(zip(LC_KEYS, LC_TYPES))


# Headers for output files
HEADERS = [TR_KEYS,
           LC_KEYS,
           ['sum'],
           FOREST_KEYS,
           ['sum'],
           ]


def read_gridlist(fname, start=0, nlines=-1, with_desc=False):

    with open(fname) as f:
        gridlist = []
        for iline, line in enumerate(f):
            if iline < start:
                continue
            if iline == start+nlines:
                break
            bits = line.split()
            lon, lat, desc = float(bits[0]), float(bits[1]), ' '.join(bits[2:])
            if with_desc:
                gridlist.append({'lon':lon, 'lat':lat, 'desc':desc})
            else:
                gridlist.append({'lon':lon, 'lat':lat})
        
    return gridlist


def get_dataset(fname, varname):

    data = Dataset(fname)
    lon = data['longitude'][:]
    lat = data['latitude'][:]
    time = data['time'][:]
    dlon = lon[1] - lon[0]
    dlat = lat[1] - lat[0]

    return data[varname], lon, lat, time, dlon, dlat


def get_site_data(data_handle, lon, lat, dlon, dlat, site):

    lon_slice = np.where( (lon >= site['lon'] - DLON_LPJG/2. + dlon/2.)
                        & (lon <= site['lon'] + DLON_LPJG/2. - dlon/2.))[0]
    lat_slice = np.where( (lat >= site['lat'] - DLAT_LPJG/2. + dlat/2.)
                        & (lat <= site['lat'] + DLAT_LPJG/2. - dlat/2.))[0]
    data_site = (data_handle[:, lat_slice, lon_slice].data).astype('uint32')

    return data_site


def make_tables(gridlist, lc_data, lon, lat, time, dlon, dlat):

    transitions_table = [] # Aggregated transitions data
    netfrac_table = [] # Net land cover fractions data
    netfrac_check_table = [] # Check that all the land cover fractions sum up to 1
    forestfrac_table = [] # Forest fractions data
    forestfrac_check_table = [] # Check that all the forest fractions sum up to 1

    for site in gridlist:
        # Retrieve slice of data for this site
        lc_data_site = get_site_data(lc_data, lon, lat, dlon, dlat, site)
        lc_previous = lc_data_site[0]
        npixels = lc_previous.size

        # Primary cells boolean matrix
        primary_cells = np.isin(lc_previous, NATURAL)
        primary_cells_left = primary_cells.any()

        # Aggregate land cover types and transition types across cell
        for year, lc_current, in zip(time[1:], lc_data_site[1:]):

            # Net fractions
            # -------------

            # Count land cover types
            row_netfrac = [site['lon'], site['lat'], year, ]
            for lc_type, lc_list in LC_TYPES.items():
                n_occurrences = np.count_nonzero(np.isin(lc_current, lc_list))
                row_netfrac.append(n_occurrences/npixels)
                # Save total number of managed forests to calc net forest fractions later
                if lc_type == 'FOREST':
                    nforests = n_occurrences
                    
            netfrac_table.append(row_netfrac)

            if NETFRAC_CHECK:
                netfrac_check_table.append(
                        [site['lon'], site['lat'], year, sum(row_netfrac[3:]), ] )

            # Transitions
            # -----------

            # Construct transition codes
            lc_previous[primary_cells] = lc_previous[primary_cells]*10 + 1 
            transitions = np.where(lc_current//100>0, lc_previous*1000, lc_previous*100) + lc_current

            # Count transition types
            row_transitions = [site['lon'], site['lat'], year, ]
            for _, tr_list in TRANSITIONS.items():
                row_transitions.append(np.count_nonzero(np.isin(transitions, tr_list))/npixels)
            transitions_table.append(row_transitions)

            # Managed forest fractions
            # ------------------------
            row_forestfrac = [site['lon'], site['lat'], year, ]
            if nforests == 0:
                row_forestfrac.extend([0., ]*NFOREST_TYPES)
            else:
                for _, mf_list in FOREST_TYPES.items():
                    row_forestfrac.append(np.count_nonzero(np.isin(lc_current, mf_list))/nforests)
            forestfrac_table.append(row_forestfrac)

            if FORESTFRAC_CHECK:
                forestfrac_check_table.append(
                        [site['lon'], site['lat'], year, sum(row_forestfrac[3:]), ] )

            # Update primary cells mask
            if primary_cells_left:
                primary_cells = np.isin(transitions, STAYS_PRIMARY)
                primary_cells_left = primary_cells.any()

            lc_previous = lc_current

    return transitions_table, netfrac_table, netfrac_check_table, \
                              forestfrac_table, forestfrac_check_table


def check_and_parse():

    usage_error = 'Usage: python hildap_tables.py [start nlines [suffix]]'
    version_error = 'Python version must be 3.7 or older.'

    try:
        assert version_info >= (3,7)
    except AssertionError:
        exit(version_error)
    
    suffix = ''
    if len(argv) == 1:
        start = nlines = None
    elif len(argv) == 3:
        try:
            start = int(argv[1])
            nlines = int(argv[2])
        except ValueError:
            exit(usage_error)
    elif len(argv) == 4:
        try:
            start = int(argv[1])
            nlines = int(argv[2])
        except ValueError:
            exit(usage_error)
        suffix = argv[3]
    else:
        exit(usage_error)

    return start, nlines, suffix


def print_info(start, nlines, fname_transitions, fname_netfrac, fname_netfrac_check,
                                                 fname_forestfrac, fname_forestfrac_check):

    if start is None:
        lines = [f"Reading all entries in gridlist.", ]
    else:
        lines = [f"Reading gridlist range [{start}:{start+nlines}].", ]

    lines.extend([ f"Transitions table written to {fname_transitions}.",
                   f"Net fractions table written to {fname_netfrac}.",
                   f"Forest fractions table written to {fname_forestfrac}.", ] )
    if NETFRAC_CHECK:
        lines.extend([ f"Net fractions check written to {fname_netfrac_check}.", ])
    if FORESTFRAC_CHECK:
        lines.extend([ f"Forest fractions check written to {fname_forestfrac_check}.", ])

    print('\n'.join(lines))


def write_table(table, fname, data_columns):

    width_coord = 10
    width_year = 6
    width_data = PRECISION + 5
    ndata_columns = len(data_columns)
    widths = [width_coord, width_coord, width_year, ] + [width_data]*ndata_columns
    format_list = [f'%{width_coord}.3f', f'%{width_coord}.3f', f'%{width_year}d', ] \
                + [f'%{width_data}.{PRECISION}f']*ndata_columns
    columns = ['lon', 'lat', 'year', ] + data_columns
    header = ''
    for col, width in zip(columns, widths):
        header += col.rjust(width)
    table = np.array(table)
    table[:,3:].round(PRECISION)
    np.savetxt(fname, table, fmt=format_list, delimiter='', header=header, comments='')


def write_output(tables, fnames, headers):

    for table, fname, header in zip(tables, fnames, headers):
        if table == []:
            continue
        print(f"Writing to {fname}...")
        write_table(table, fname, header)


def main(fname_dataset, fname_gridlist, start=None, nlines=None, suffix=''):

    suffix = f'.{suffix}' if suffix != '' else ''
    fname_transitions = 'hildaplus_gross_1901_2020.txt' + suffix
    fname_netfrac = 'hildaplus_netfrac_1901_2020.txt' + suffix
    fname_netfrac_check = 'hildaplus_netfrac_check_1901_2020.txt' + suffix
    fname_forestfrac = 'hildaplus_forestfrac_1901_2020.txt' + suffix
    fname_forestfrac_check = 'hildaplus_forestfrac_check_1901_2020.txt' + suffix
    fnames = fname_transitions,                         \
             fname_netfrac, fname_netfrac_check,        \
             fname_forestfrac, fname_forestfrac_check

    # Print runtime information to the screen
    print_info(start, nlines, *fnames)

    # Read LPJ-GUESS gridlist
    if start is None or nlines is None:
        gridlist = read_gridlist(fname_gridlist)
    else:
        gridlist = read_gridlist(fname_gridlist, start=start, nlines=nlines)

    # Retrieve dataset
    print("Retrieving data...")
    data, lon, lat, time, dlon, dlat = get_dataset(fname_dataset, 'LULC_states')

    # Construct tables
    print("Constructing tables...")
    tables = make_tables(gridlist, data, lon, lat, time, dlon, dlat)

    # Write output
    write_output(tables, fnames, HEADERS)
    

if __name__ == '__main__':

    # Check system requirements and args correctness
    start, nlines, suffix = check_and_parse()

    # Start timer
    t_start = timer()

    # Run main function
    main(FNAME_DATAFILE, FNAME_GRIDLIST, start, nlines, suffix)

    # Stop timer
    t_end = timer()
    print('Elapsed time (s): ', round(t_end - t_start, 3))

