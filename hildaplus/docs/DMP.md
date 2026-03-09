# hildaplus (formerly hildaplus_global_05deg_lpjguess)

# Authors

Martin Wittenbrink
(martin.wittenbrink@kit.edu) 
ORCID: 0000-0002-6227-7609

David Martín Belda
(david.belda@kit.edu)
ORCID: 0009-0005-7702-8233


# Project description

LPJ-GUESS landcover input files derived from HILDA+. The original dataset has been smoothed to eliminate a flickering artifact (very fast land cover changes) and upscaled to a 0.5 degree resolution.

# Data   
## Data re-used 

- HILDA+ Landuse/Landcover data set, updated version that is not yet publicly available. We used the LULC_states and the forest_management files.
- location: /pd/data/lpj/input_data/Hilda/Hildaplus_vGLOB-2.0/hildap_vGLOB-2.0-forest_netcdf/
- The LULC_states and forest management files were merged into one dataset and then smoothed because of flickering issues.
- The dataset was subsequently up-scaled to 0.5deg resolution.

## Data produced

Landcover net fractions, gross transitions and lu_forest input files for LPJ-GUESS on 0.5 degree spatial resolution were created based on HILDA+ (text format)


- `data/hildaplus_gross_1901_2020.txt`: Gross transitions between land covers following Anita's [work](https://esd.copernicus.org/articles/8/91/2017/) with Hilda.
- `data/hildaplus_netfrac_1901_2020.txt`: Net land cover fractions following Anita's work.
- `data/hildaplus_forestfrac_1901_2020.txt`: Managed forest cover fractions 

## Technical classification

- Expected size of the dataset: Less than 10 GB
- Format: text files

## Which tools, software, or processes are used to generate or collect data? 

- Python 3 with:
    - Numpy
    - netCDF4

## Version Control

Remote repository: 
https://gitlab.imk-ifu.kit.edu/lemg/wittenbrink-m/hildaplus_global_05deg_lpjguess

# Requirements 

The requirements needed to run the project are defined in hildaplus/requirements.txt

# Long-term preservation  

The project data should be stored for at least 5 years after the end of the project. The project data consist of final runs (outputs) used in publications/reports, input used to run the simulations, and in some cases, some intermediate runs.

##  Storage location

/bg/data/lpj/LPJ-GUESS/input/LU
       
# Backup

The backups will be made every month. The responsible person is Martin Wittenbrink.

# Data sharing and re-use

The data can be shared internally with LEMG members,
as long as they don't pass on the data.
The data might be shared externally only with 
individual approval.

