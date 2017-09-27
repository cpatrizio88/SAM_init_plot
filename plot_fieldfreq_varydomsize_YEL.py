import matplotlib as mpl
mpl.use('Agg')
from netCDF4 import Dataset
import site
site.addsitedir('/glade/scratch/patrizio/thermolib/')
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
#import SAM_init_plot.block_fns
#import SAM_init_plot.misc_fns
#import thermolib
#from SAM_init_plot.block_fns import blockave2D, blockave3D, blockxysort2D
from block_fns import blockave2D,blockave3D
#from thermolib.constants import constants
from constants import constants
from wsat import wsat
#from thermolib.wsat import wsat

matplotlib.rcParams.update({'font.size': 26})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 2})
matplotlib.rcParams.update({'legend.fontsize': 22})


plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to170_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday90to160_64vert_ubarzero_FREQDIST/'

fout = '/Users/cpatrizio/Google Drive/MS/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_FREQDIST/'

fpath2D = '/glade/scratch/patrizio/OUT_2D_nc/'
fpath3D = '/glade/scratch/patrizio/OUT_3D_nc/'
fout = '/glade/scratch/patrizio/OUT_2D_FIGS/'
fpathSTAT = '/glade/scratch/patrizio/OUT_STAT_nc/'

#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

nc_in2D1 = glob.glob(fpath2D + '*256x256*3000m*day230to250*302K.nc')[0]
nc_in3D1 = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]

nc_in2D2 = glob.glob(fpath2D +  '*512x512*3000m*day180to195*302K.nc')[0]
nc_in3D2 = glob.glob(fpath3D +  '*512x512*3000m*day180to195*302K.nc')[0]

nc_in2D3 = glob.glob(fpath2D + '*1024x1024*3000m*day170to180*302K.nc')[0]
nc_in3D3 = glob.glob(fpath3D +  '*1024x1024*3000m*day170to180*302K.nc')[0]

nc_in2D4 = glob.glob(fpath2D + '*2048*3000m*.nc')[0]
nc_in3D4 = glob.glob(fpath3D + '*2048*3000m*.nc')[0]

nc_fs = [nc_in2D1, nc_in2D2, nc_in2D3, nc_in2D4]
nc_3Dfs = [nc_in3D1, nc_in3D2, nc_in3D3, nc_in3D4]

#domsize=768
#domsize=1536
#domsize=3072

domsizes = [768, 1536, 3072, 6144]
colors=['k', 'r', 'g', 'm']

#averaging time period for 2D & 3D fields
ntave2D=24
ntave3D=4

print 'loading netCDF files'

fields = np.array([])

db=1

varname='W'

t=-1
nave=5
aveperiod = 24*nave

for i, nc_in2D in enumerate(nc_fs):
#nc_data3D = Dataset(nc_in3D) 
    nc_in3D = nc_3Dfs[i]
    print 'loading', nc_in2D
    print 'loading', nc_in3D
    domsize = domsizes[i]
    if domsize == 768:
        nave=8
    elif domsize == 1536:
        nave=8
    else:
        nave=8
    aveperiod = 24*nave
    aveperiod3D = 4*nave
    nc_data2D = Dataset(nc_in2D)
    nc_data3D = Dataset(nc_in3D)
    #varis3D = nc_data3D.variables
    varis2D = nc_data2D.variables
    varis3D = nc_data3D.variables
    x = varis2D['x'][:]
    y = varis2D['y'][:]
    #z = varis3D['z'][:]
    #p = varis3D['p'][:]
    #p = p*1e2
    
    #nt3D = t3D.size
    #nz = z.size
    nx = x.size
    ny = y.size
    if varname == 'W':
        units = varis3D[varname].units.strip()
        vari = varis3D[varname][t-aveperiod3D:t,:,:,:]
        #vari = np.mean(np.mean(vari,axis=0),axis=0)
    else:
        units = varis2D[varname].units.strip()
        vari = varis2D[varname][t-aveperiod:t,:,:]
        vari = np.mean(vari,axis=0)
    
    t3D = varis2D['time'][:]
    t2D = varis2D['time'][:]
    nt2D = t2D.size
    #vari = varis2D[varname]
    #field = np.mean(vari[t-aveperiod:t,:,:], axis=0)
    #field = blockave2D(vari, db)
    field=vari
    nxprime = nx
    nyprime = ny
    nz = field.shape[1]
    #nxprime = nx / db
    #nyprime = ny / db
    nt = field.shape[0]
    field = field.flatten()
    bins = np.linspace(-25,25, 50)
    #nt = field.shape[
    #bins=50
    
    N = nz*nxprime*nyprime*nt
    
    #freqs, bin_edges = np.histogram(field, bins=bins, weights=np.zeros_like(field) + 1. / (nz*nxprime*nyprime*nt))
    #freqs, bin_edges = np.histogram(field, bins=bins)
    #bin_c = (bin_edges[1:] + bin_edges[:-1])/2
    
    print 'variance of {:s} = {:3.4f}'.format(varname, np.var(field))
    print 'max of {:s} = {:3.4f}'.format(varname, np.max(field))
        
    plt.figure(1)
    plt.hist(field, bins=bins, weights=np.zeros_like(field) + 1. / (nz*nxprime*nyprime*nt), histtype='step', lw=2.5, color=colors[i], label='{:d} km, day {:3.0f} to {:3.0f}, N = {:2.2e}'.format(domsizes[i], t2D[t-aveperiod], t2D[t], N))
    #plt.plot(bin_c, freqs, 'x-', color=colors[i], label='{:d} km, day {:3.0f} to {:3.0f}, N = {:2.2e}'.format(domsizes[i], t2D[t-aveperiod], t2D[t], N))
    plt.yscale('log', nonposy='clip')
    plt.savefig(fout + '{:s}freqNEW_db{:d}.pdf'.format(varname, db))

if varname == 'W':
    #titlename = 'vertically averaged W'
    titlename = 'W'
else:
    titlename = varname
    
if db == 1:
    plt.title(r'$w$ frequency distribution'.format(titlename))
else:
    plt.title(r'$w$ frequency distribution, block-averaging over ({:2.0f} km)$^2$'.format(titlename, db*(np.diff(x)[0])/1e3))
    
if (varname == 'ZC') or (varname == 'ZE'):
   plt.ylim(0,0.03)
if varname == 'W500':
   plt.xlim(0, 1)
   plt.ylim(0, 0.1)

plt.xlabel('{:s} ({:s})'.format(varname, units))
plt.ylabel('frequency')
#plt.ylim((0,0.2))
#plt.xlim((0,0.02))
plt.axvline(0, color='k', alpha=0.5, linewidth=0.5)
plt.xlim((-25,25))
plt.legend(loc='best')
plt.savefig(fout + '{:s}freqNEW_db{:d}.pdf'.format(varname, db))
plt.close()




    
    