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
#from SAM_init_plot.misc_fns import radprof3D, radprof
#from misc_fns import radprof3D, radprof
#from SAM_init_plot.block_fns import blockave2D, blockave3D, blockxysort2D
from block_fns import blockave2D, blockave3D
from constants import constants
import matplotlib.colors as colors

matplotlib.rcParams.update({'font.size': 28})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 2})
matplotlib.rcParams.update({'legend.fontsize': 24})

plt.style.use('seaborn-white')

#fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
#fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to170_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/MS/figures/SST302/3072km_SAM_aggrday90to160_64vert_ubarzero_FREQDIST/'

fpath2D = '/glade/scratch/patrizio/OUT_2D_nc/'
fpath3D = '/glade/scratch/patrizio/OUT_3D_nc/'
fout = '/glade/scratch/patrizio/FREQDIST_FIGS/'

#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*512x512*3000m*day090to140*302K.nc')[0]

#nc_in2D = glob.glob(fpath2D + '*1024x1024*3000m*day090to130*302K.nc')[0]

#nc_fs = glob.glob(fpath2D + '*256*256*3000m*302K.nc')
#nc_fs = glob.glob(fpath2D + '*512*512*3000m*302K.nc')
#nc_fs = glob.glob(fpath2D + '*1024x1024*3000m*302K.nc')
nc_fs = nc_in2D = glob.glob(fpath2D + '*2048*3000m*.nc')[0]
nc_fs = nc_fs[:1]
#domsize=768
#domsize=1536
#domsize=3072
domsize=6144

nc_data2D = Dataset(nc_in2D)
varis2D = nc_data2D.variables

x = varis2D['x'][:]
y = varis2D['y'][:]
#z = varis3D['z'][:]
#p = varis3D['p'][:]
#p = p*1e2

tstart = 0 

#nt3D = t3D.size
#nz = z.size
nx = x.size
ny = y.size

#averaging time period for 2D & 3D fields
ntave2D=24
ntave3D=4

print 'loading netCDF files'

fields = np.array([])
times = np.array([])
db=1

tinit = 200

varname='PW'

for nc_in2D in nc_fs:
#nc_data3D = Dataset(nc_in3D) 
    print 'loading', nc_in2D
    nc_data2D = Dataset(nc_in2D)
    #varis3D = nc_data3D.variables
    varis2D = nc_data2D.variables
    t2D = varis2D['time'][:]
    nt2D = t2D.size
    field = varis2D[varname][:]
    #ntrunc = nt2D%ntave2D
    #ntrunc=0
    #field = field[ntrunc:,:,:]
    #print field.shape
    units = varis2D[varname].units
    #field_tmp = field.reshape(nt2D/ntave2D, ntave2D, nx, ny)
    #field_tave = np.mean(field_tmp, axis=1)
    if db == 1:
        #fields_temp = field_tave
        fields_temp = field
    else:
        fields_temp = blockave3D(field, db)
    #        buoyflux = np.vstack((buoyflux, buoyflux_temp)) if buoyflux.size else buoyflux_temp
    fields = np.vstack((fields, fields_temp)) if fields.size else fields_temp
    #times_temp = np.arange(
    #times = np.concatenate((times, t2D)) if times.size else t2D

ntimes = fields.shape[0]

times = np.arange(tinit, tinit+ntimes)

nx = fields.shape[1]

#nfieldbins= 1024
#nfieldbins = np.linspace(-0.1, 0.1, 300)

nfieldbins=1024

#nfieldbins = 1000

#nfieldbins=1024

fields = fields.reshape(ntimes, nx*nx)

print 'contouring'
 
#calculate 2D (CRH rank, time) map of CRH frequency  ###
#times = np.arange(tstart, tstart+ntimes)

fieldss, tt =  np.meshgrid(np.arange(nx*nx), times)

tcoords = tt.flatten()
fieldcoords = fields.flatten()
fieldfreqs, xedges, yedges = np.histogram2d(np.transpose(tcoords), fieldcoords, bins=[ntimes, nfieldbins], weights=np.zeros_like(fieldcoords) + 1. / (nx*nx))
#fieldfreqs, xedges, yedges = np.histogram2d(np.transpose(tcoords), fieldcoords, bins=[ntimes, nfieldbins])

xmids = (xedges[1:] + xedges[:-1])/2.
ymids = (yedges[1:] + yedges[:-1])/2.

xxmids, yymids = np.meshgrid(xmids, ymids)

#extent = [-0.1, 0.3, yedges[-1], yedges[0]]

#vals = np.arange(0, 200, 5)

#ymin=0
#ymax=0.5*np.round(np.max(fields))

#ymin=-0.1
#ymax=0.1

plt.figure(1)
#plt.pcolormesh(xxmids, yymids, np.transpose(fieldfreqs), norm=colors.LogNorm(vmin=1, vmax=fieldfreqs.max()), cmap=cm.inferno)
#plt.pcolormesh(xxmids, yymids, np.transpose(fieldfreqs), vmin=0, vmax=0.04, cmap=cm.inferno)
plt.pcolormesh(xxmids, yymids, np.transpose(fieldfreqs), cmap=cm.BuPu)
#plt.ylim(ymin, ymax)
cb = plt.colorbar()
cb.set_label('frequency')

if db == 1:
    tt1=plt.title(r'daily-averaged {:s} frequency distribution, domain size = ({:3.0f} km)$^2$'.format(varname, domsize))
else:
    tt1 = plt.title(r'daily-averaged {:s} frequency distribution, domain size = ({:3.0f} km)$^2$, block-averaging over ({:2.0f} km)$^2$'.format(varname, domsize, db*(np.diff(x)[0])/1e3))
tt1.set_position([0.5, 1.008]) 
plt.ylabel('{:s} ({:s})'.format(varname, units))
plt.xlabel('time (days)')
plt.savefig(fout + '{:s}freqvst_db{:d}new.png'.format(varname, db))
plt.close()




    
    