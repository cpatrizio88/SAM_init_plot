from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from SAM_init_plot.block_fns import blockave2D, blockave3D
from thermolib.constants import constants
import matplotlib.colors as colors

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fpathSTAT =  '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to170_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday90to160_64vert_ubarzero_FREQDIST/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STAT/'

#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

#nc_in2D = glob.glob(fpath2D + '*1024x1024*3000m*day090to130*302K.nc')[0]

nc_f2D1 = glob.glob(fpath2D + '*256x256*3000m*302K.nc')
nc_f2D2 = glob.glob(fpath2D + '*512x512*3000m*302K.nc')
nc_f2D3 = glob.glob(fpath2D + '*1024x1024*3000m*302K.nc')

nc_fSTAT1 = glob.glob(fpathSTAT + '*256x256*3000m*250days_302K.nc')[0]
nc_fSTAT2 = glob.glob(fpathSTAT + '*512x512*3000m*195days_302K.nc')[0]
nc_fSTAT3 = glob.glob(fpathSTAT + '*1024x1024*3000m*200days_302K.nc')[0]

domsize=768
#domsize=1536
#domsize=3072

nc_f2Ds = [nc_f2D1, nc_f2D2, nc_f2D3]
nc_fSTATs = [nc_fSTAT1, nc_fSTAT2, nc_fSTAT3]
domsizes=[768, 1536, 3072]
colors = ['k', 'r', 'g']

tstarts = [0,90,90]

for i, nc_f in enumerate(nc_fSTATs):
        
    print 'domsize', domsizes[i]
    
    vals = np.array([])
    
    times = np.array([])
    
    varname='QV'


    print 'loading', nc_f
    
    #averaging time period for 2D & 3D fields
    ntave2D=24
    ntave3D=4
    nc_dataSTAT = Dataset(nc_f)
    varisSTAT = nc_dataSTAT.variables
    t2D = varisSTAT['time'][:]
    z = varisSTAT['z'][:]
    nt2D = t2D.size
    vari = varisSTAT[varname]
    field = vari[:]
    
    z_BLi = np.where(z > 1.5e3)[0][0]
    
    
    field_BLave = np.mean(field[:,:z_BLi], axis=1)
    

    units = vari.units
    
    vals = np.concatenate((vals, field_BLave)) if vals.size else field_BLave
    times = np.concatenate((times, t2D)) if times.size else t2D
        
    print 'plotting'
    w = 24*5
    
    vals_smooth = moving_average(vals, w)
    
    plt.figure(1)
    plt.plot(times, vals, color=colors[i], alpha=0.5, label='{:d} km'.format(domsizes[i]))
    plt.plot(times[w/2:-w/2+1], vals_smooth, color=colors[i], linewidth=1.5)
    plt.xlabel('time (days)')
    plt.ylabel('{:s} ({:s})'.format(varname, vari.units.strip()))
    plt.title('average {:s} between z = 0 and z ={:3.1f} km'.format(varname, z[z_BLi]/1e3))
    plt.savefig(fout + 'BL{:s}250days.pdf'.format(varname))
  

plt.legend()
plt.savefig(fout + 'BL{:s}250days.pdf'.format(varname))
plt.close()
        

    