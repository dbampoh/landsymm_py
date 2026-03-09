# HILDA+ UPSCALING

These scripts use the pre-processed, smoothed Hilda+ Land Cover file created by Martin (martin.wittenbrink@kit.edu) to produce upscaled versions usable by LPJ-GUESS. Python 3.7 or higher is required!

## Files:

* `hildap_tables.py`: this one creates tables of land cover net fractions, transitions and managed forest relative cover, and saves them to text files.
* `launch_tables.py`: this one is to run the above script in the IFU cluster (keal)
* `README.md`: this file


## Usage:

### Running in serial (e.g. your laptop)

1. Edit the `hildap_tables.py`
    - `FNAME_DATAFILE`: path/fname to the smoothed file
    - `DLON_LPJG`: Upscaled longitude resolution
    - `DLAT_LPJG`: Upscaled latitude resolution
    - `PRECISION`: Output is rounded to `PRECISION` decimals
    - `FNAME_GRIDLIST`: Path/name of the gridlist file
    - `NETFRAC_CHECK`: If True, it creates a file adding up the net fractions (they all should add up to 1)
    - `FORESTFRAC_CHECK`: If True, it creates a file adding up the managed forest fractions (if there are managed forests in the gridcell, they all should add up to 1. Otherwise 0.)

2. To run, enter `python hildap_tables.py [start nlines [suffix]]`
    - If executed without arguments, the script will read the full gridlist
    - `start`: [integer] starting gridlist line (0 = first line)
    - `nlines`: [integer] number of lines to read
    - `suffix`: [string] adds ".suffix" to the output files

3. Output:
    - `hildaplus_gross_1901_2020.txt`: Gross transitions between land covers following Anita's [work](https://esd.copernicus.org/articles/8/91/2017/) with Hilda.
    - `hildaplus_netfrac_1901_2020.txt`: Net land cover fractions following Anita's work.
    - `hildaplus_netfrac_check_1901_2020.txt`: Sum of net fractions for every gridcell/year. This should always add up to 1 within machine precision limits.
    - `hildaplus_forestfrac_1901_2020.txt`: Managed forest cover fractions 
    - `hildaplus_forestfrac_check_1901_2020.txt`: Sum of forest fractions for every gridcell/year. This should always add up to 1, within machine precision limits, or 0, when there are no managed forests in the gridcell.

4. Examples:
    - `python hildap_tables.py`: Processes the whole gridlist
    - `python hildap_tables.py 3 20`: Processes 20 lines of the gridlist, starting at line 4
    - `python hildap_tables.py 3 20 test`: Processes 20 lines of the gridlist, starting at line 4, and adds the extension `.test` to the output files


### Running on keal

1. Edit the `launch_tables.py` script:
    - `NODES`: Nodes to use to run the script. You should introduce only the node numbers in the form of a python list. All the nodes should be in the same partition.
    - `NODE_IO`: Node to submit the post-processing job.
    - `WT_DAYS[HOURS, MINUTES, SECONDS]`: Time limit
    - `JOBNAME`: Job name
    - `FNAME_GRIDLIST`: Gridlist file path/name
    - `TEST`: If True, the submit and append (postprocessing) scripts will be generated, but not submitted to the cluster.

2. Run `python launch_tables.py`

3. If everything went well, you should have the output files mentioned above, plus a directory 'outfiles' with the un-appended files. This directory can be removed when the job is finished.

## Technical information 

LPJ-GUESS classifies LCTs as:

- b: Barren
- c: Cropland
- p: Pasture
- v: Natural (primary)
- s: Natural (secondary)
- f: Managed forest

Hilda+v2.0 classifies LCTs as:

-  0: Ocean
- 11: Urban
- 22: Cropland
- 33: Pasture
- 40-45: Forest (natural or managed)
- 55: Grassland/Shrubland
- 66: Other land (barren, etc...)
- 77: Water
- 99: No data

A transition between LCTs is represented by concatenating the *previous* and the *current* LCT codes. For example, 3322 represents a transition from pasture to cropland in Hilda+. A code 'vc' represents a primary natural land converted to cropland in LPJ-GUESS.

### Land cover type equivalences

Martin's smoothed dataset uses the usual Hilda+ land cover codes, and includes the forest management information in the same file by multiplying the codes of managed forests by 10. I will classify ocean, water and no data gridcells as barren. The resuling equivalence table reads:

| Code    | LPJ-GUESS LCT                           |
|---------|-----------------------------------------|
| 0       | Barren (ocean)                          |
| 11      | Urban                                   |
| 22      | Crop                                    |
| 33      | Pasture                                 |
| 40-45   | Natural (unmanaged forest)              |
| 400-450 | Forest (managed forest)                 |
| 55      | Natural (grass/shrubland)               |
| 66      | Barren (barren/other land)              |
| 77      | Barren (water)                          |
| 99      | Barren (no data)                        |

Note that the LPJ-GUESS gridlist contains only land points (or it should). So counting ocean as barren should not make a big difference, as it is within the error margin of the target resolution.

### Land cover transitions

To characterize the transitions, the following symbols for Hilda+ LCTs are introduced:

- X, Y: Any LCT except for Forest or Grassland/Shrubland (i.e., one of [11, 22, 33, 66]).
- F: Forest. One of [40, 41, 42, 43, 44, 45].
- G: Grassland/Shrubland: [55].

Types 0, 77 and 99 do not change across the dataset, so they are not considered here. Unmanaged forests and grasslands present at the beginning of the dataset are considered primary. If the previous LCT is primary, A flag `1` is appended to its code before concatenating it with the current LCT code.

Examples:

- A non-primary natural forest with code 44 that switches to managed 440 would have transition code 44440.

- The transition code for a primary grassland that switches to unmanaged forest of type 43 is 55143.

The following table summarizes the equivalences between all possible transition codes so constructed and their LPJ-GUESS gross transition category:

| Code | Transition |
|------|------------|
| YX   | yx         |
| XF   | xs         |
| XF0  | xf         |
| XG   | xs         |
| FX   | sx         |
| F0X  | fx         |
| F1X  | vx         |
| FF   | -          |
| F0F  | fs         |
| F1F  | -          |
| FF0  | sf         |
| F0F0 | -          |
| F1F0 | vf         |
| GF   | -          |
| GF0  | sf         |
| G1F  | -          |
| G1F0 | vf         |
| FG   | -          |
| F0G  | fs         |
| F1G  | -          |
| GG   | -          |
| G1G  | -          |
| GX   | sx         |
| G1X  | vx         |

Inverting this table yields the transformed Hilda+ codes that go into each LPJ-GUESS transition category:

| Transition | Code       |
|------------|------------|
| yx         | YX         |
| xs         | XF, XG     |
| xf         | XF0        |
| sx         | FX, GX     |
| fx         | F0X        |
| vx         | F1X, G1X   |
| sf         | FF0, GF0   |
| fs         | F0F, F0G   |
| vf         | F1F0, G1F0 |
