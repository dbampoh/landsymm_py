# hildaplus_smoothing

Set of scripts to reduce flickering issues in HILDA+ dataset with gaussian filtering.

smoothingscript_run.py contains the functions to perform the smoothing, split the data into manageable chunks, run it prarallelized and store the results in separate files

create_netcdf_from_slices.py creates a netcdf file from the smoothed data slices

adapt_nc.py sets a netcdf description afterwards
