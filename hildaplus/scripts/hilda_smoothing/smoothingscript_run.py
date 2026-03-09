import numpy as np
import netCDF4 as nc
from scipy import ndimage as nd
import joblib
from joblib import Parallel, delayed


def hildafomafusion(hildastates,hildafoma):				# function to merge the hilda+ LULC states and the forest management data into one dataset
    return hildastates*(np.ones_like(hildafoma)+(hildafoma==1)*9)     	# pixels where forest management is true are mulitplied by 10


def reduce_flickering_gaussian_1pixel(slice1,sigma=3.0):
    uniq=np.unique(slice1,return_inverse=True)    			# get different values in timeseries
    if len(uniq[0])>1:                            			# only continue, if more than 1 value (to save computing time)
        u=[uniq[1]==i for i in range(len(uniq[0]))]			# creates sepaarate timeseries for all LULC classes in pixel
        gf=[nd.gaussian_filter(ui.astype(float),sigma,mode='constant',cval=0) for ui in u]		# gaussian filtering
        imax=np.argmax(gf,axis=0)					# now construct a new timeseries maximum categories from gaussian filtering
        newslice=slice1.copy()
        for i in range(len(uniq[0])):
            newslice[imax==i]=uniq[0][i]
        return newslice
    else:
        return slice1
    
def apply_gaussian_3times(slice1,sigma=3.0):				# apparently applying it once isn't enough
    return reduce_flickering_gaussian_1pixel(reduce_flickering_gaussian_1pixel(reduce_flickering_gaussian_1pixel(slice1,sigma),sigma),sigma)

def runlonslice(hdata):							# run the smoothing for a given block of Hilda+ data
    return [[apply_gaussian_3times(hdata[:,i,j]) for i in range(hdata.shape[1])] for j in range(hdata.shape[2])]

def run_it(i_slice):							# collect the data and run the smoothing for a dataset of 121x18000x100 pixels (years x latitude x longitude, 1/360 of the whole global data)
    # i_slice is an index running from 0 to 359 in order to cut the global process into 360 parallelizable processes
    hildastatesfile='/pd/data/lpj/input_data/Hilda/Hildaplus_vGLOB-2.0/hildap_vGLOB-2.0-forest_netcdf/hildaplus_GLOB-2-0-f_states.nc'
    hildafomafile='/pd/data/lpj/input_data/Hilda/Hildaplus_vGLOB-2.0/hildap_vGLOB-2.0-forest_netcdf/hildaplus_GLOB-2-0-f_forest-management.nc'
    hstates=nc.Dataset(hildastatesfile,'r')
    hfoma=nc.Dataset(hildafomafile,'r')
    hstatesdata=hstates['LULC_states'][:,:,i_slice*100:(i_slice*100+100)].data
    hfomadata=hfoma['forest_management'][:,:,i_slice*100:(i_slice*100+100)].data
    hildadata=hildafomafusion(hstatesdata,hfomadata)
    slicedata=runlonslice(hildadata)
    np.savez_compressed('/pd/data/lpj/wittenbrink-m/smoothing/slice_npys/slice_'+str(i_slice)+'.npz', a=slicedata)  # save the smoothed data as a compressed npy file
    return

Parallel(n_jobs=20)([delayed(run_it)(i) for i in range(360)])		# run it in parallel 
                                  
