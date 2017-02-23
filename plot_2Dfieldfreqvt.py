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

#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

#nc_in2D = glob.glob(fpath2D + '*1024x1024*3000m*day090to130*302K.nc')[0]

nc_fs = glob.glob(fpath2D + '*256x256*3000m*302K.nc')
domsize=768
#domsize=1536
#domsize=3072

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

db=8

varname='CWP'

for nc_in2D in nc_fs:
#nc_data3D = Dataset(nc_in3D) 
    print 'loading', nc_in2D
    nc_data2D = Dataset(nc_in2D)
    #varis3D = nc_data3D.variables
    varis2D = nc_data2D.variables
    t2D = varis2D['time'][:]
    nt2D = t2D.size
    field = varis2D[varname][:]
    
    ntrunc = nt2D%ntave2D
    #ntrunc=0
    field = field[ntrunc:,:,:]
    units = varis2D[varname].units
    field_tmp = field.reshape(nt2D/ntave2D, ntave2D, nx, ny)
    field_tave = np.mean(field_tmp, axis=1)
    if db == 1:
        fields_temp = field_tave
    else:
        fields_temp = blockave3D(field_tave, db)
    #        buoyflux = np.vstack((buoyflux, buoyflux_temp)) if buoyflux.size else buoyflux_temp
    fields = np.vstack((fields, fields_temp)) if fields.size else fields_temp

ntimes = fields.shape[0]

nx = fields.shape[1]

nfieldbins= 1024

fields = fields.reshape(ntimes, nx*nx)

print 'contouring'
 
#calculate 2D (CRH rank, time) map of CRH frequency  ###
times = np.arange(tstart, tstart+ntimes)

fieldss, tt =  np.meshgrid(np.arange(nx*nx), times)

tcoords = tt.flatten()
fieldcoords = fields.flatten()
fieldfreqs, xedges, yedges = np.histogram2d(np.transpose(tcoords), fieldcoords, bins=(ntimes, nfieldbins), weights=np.zeros_like(fieldcoords) + 1. / (nx*nx))

xmids = (xedges[1:] + xedges[:-1])/2.
ymids = (yedges[1:] + yedges[:-1])/2.

xxmids, yymids = np.meshgrid(xmids, ymids)

#extent = [-0.1, 0.3, yedges[-1], yedges[0]]

#vals = np.arange(0, 200, 5)

ymin=0
ymax=0.5*np.round(np.max(fields))

plt.figure(1)
#plt.pcolormesh(xxmids, yymids, np.transpose(fieldfreqs), norm=colors.LogNorm(vmin=1, vmax=fieldfreqs.max()), cmap=cm.inferno)
plt.pcolormesh(xxmids, yymids, np.transpose(fieldfreqs), vmin=0, vmax=0.05, cmap=cm.inferno)
plt.ylim(ymin, ymax)
cb = plt.colorbar()
cb.set_label('frequency')

if db == 1:
    plt.title(r'{:s} frequency distribution, domain size = ({:3.0f} km)$^2$'.format(varname, domsize))
else:
    plt.title(r'{:s} frequency distribution, domain size = ({:3.0f} km)$^2$, block-averaging over ({:2.0f} km)$^2$'.format(varname, domsize, db*(np.diff(x)[0])/1e3))
plt.ylabel('{:s} ({:s})'.format(varname, units))
plt.xlabel('time (days)')
plt.savefig(fout + '{:s}freqvst_db{:d}new.pdf'.format(varname, db))
plt.close()




    
    