from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from SAM_init_plot.block_fns import blockave2D, blockave3D
from thermolib.constants import constants
import matplotlib.colors as colors

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to170_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday90to160_64vert_ubarzero_FREQDIST/'

fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_FREQDIST/'

#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

nc_in2D1 = glob.glob(fpath2D + '*256x256*3000m*day230to250*302K.nc')[0]

nc_in2D2 = glob.glob(fpath2D +  '*512x512*3000m*day180to195*302K.nc')[0]

nc_in2D3 = glob.glob(fpath2D + '*1024x1024*3000m*day180to190*302K.nc')[0]

nc_fs = [nc_in2D1, nc_in2D2, nc_in2D3]

#domsize=768
#domsize=1536
#domsize=3072

domsizes = [768, 1536, 3072]
colors=['k', 'r', 'g']

#averaging time period for 2D & 3D fields
ntave2D=24
ntave3D=4

print 'loading netCDF files'

fields = np.array([])

db=1

varname='W500'

t=-1
nave=5
aveperiod = 24*nave

for i, nc_in2D in enumerate(nc_fs):
#nc_data3D = Dataset(nc_in3D) 
    print 'loading', nc_in2D
    nc_data2D = Dataset(nc_in2D)
    #varis3D = nc_data3D.variables
    varis2D = nc_data2D.variables
    x = varis2D['x'][:]
    y = varis2D['y'][:]
    #z = varis3D['z'][:]
    #p = varis3D['p'][:]
    #p = p*1e2
    
    #nt3D = t3D.size
    #nz = z.size
    nx = x.size
    ny = y.size
    t2D = varis2D['time'][:]
    nt2D = t2D.size
    vari = varis2D[varname]
    field = np.mean(vari[t-aveperiod:t,:,:], axis=0)
    field = blockave2D(field, db)
    nxprime = nx / db
    nyprime = ny / db
    
    field = field.flatten()
    
    freqs, bin_edges = np.histogram(field, bins=50, weights=np.zeros_like(field) + 1. / (nxprime*nyprime))
    bin_c = (bin_edges[1:] + bin_edges[:-1])/2
        
    plt.figure(1)
    plt.plot(bin_c, freqs, 'x-', color=colors[i], label='{:d} km, day {:3.0f} to {:3.0f}'.format(domsizes[i], t2D[t-aveperiod], t2D[t]))

    
if db == 1:
    plt.title(r'{:s} frequency distribution'.format(varname))
else:
    plt.title(r'{:s} frequency distribution, block-averaging over ({:2.0f} km)$^2$'.format(varname, db*(np.diff(x)[0])/1e3))
    
if (varname == 'ZC') or (varname == 'ZE'):
   plt.ylim(0,0.03)
if varname == 'W500':
   plt.xlim(0, 1)
   plt.ylim(0, 0.1)

plt.xlabel('{:s} ({:s})'.format(varname, vari.units.strip()))
plt.ylabel('frequency')
plt.legend()
plt.savefig(fout + '{:s}freq_db{:d}.pdf'.format(varname, db))
plt.close()




    
    