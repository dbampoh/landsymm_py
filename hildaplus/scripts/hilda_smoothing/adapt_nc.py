import numpy as np
import netCDF4 as nc

### the purpose of this script is to add a description to the smoothed HILDA+ netcdf file ###


filename_oldhilda='/pd/data/lpj/input_data/Hilda/Hildaplus_vGLOB-2.0/hildap_vGLOB-2.0-forest_netcdf/hildaplus_GLOB-2-0-f_states.nc'
filename='hildap_smoothed3.nc'

ds4=nc.Dataset(filename_oldhilda)

ncf=nc.Dataset(filename,'r+')
time=ncf.createVariable('time','uint16',('time'))
lat=ncf.createVariable('latitude','float32',('lat'))
lon=ncf.createVariable('longitude','float32',('lon'))
time[:]=ds4['time'][:].data
lat[:]=ds4['latitude'][:].data
lon[:]=ds4['longitude'][:].data
ncf.description='HILDA+ land use/cover reconstruction - LULC states - SMOOTHED'
ncf.producer='Martin Wittenbrink'
ncf.close()
ds4.close()
