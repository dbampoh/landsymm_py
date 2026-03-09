import numpy as np
import netCDF4 as nc

### Script that constructs a netcdf file from smaller npy-files of smoothed hilda+ data ###


slicepath='/pd/data/lpj/wittenbrink-m/smoothing/slice_npys/'		# path were the smaller npy-files to construct the netcdf are stored
ncfilename='hildap_smoothed2.nc'					# filename of the resulting netcdf file

f = nc.Dataset(ncfilename,'w')
f.createDimension('time',121)
f.createDimension('lat',18000)
f.createDimension('lon',None)

newstates=f.createVariable('LULC_states','uint16',('time','lat','lon'),zlib=True)
for i in range(360):
    newvals=f['LULC_states']
    print(i)
    sli=np.load(f'{slicepath}slice_{i}.npz')['a']
    for j in range(100):
        newvals[:,:,i*100+j]=sli.swapaxes(0,2)[:,:,j]
    del sli
    del newvals
f.close()
