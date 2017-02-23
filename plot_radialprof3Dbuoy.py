from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
import SAM_init_plot.block_fns
import SAM_init_plot.misc_fns
from SAM_init_plot.misc_fns import radprof3D
from SAM_init_plot.block_fns import blockave2D, blockave3D, blockxysort2D
from thermolib.constants import constants

c=constants()

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fin = '/Users/cpatrizio/data/SST302/768km/'
fout = '/Users/cpatrizio/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_RADIALXSECTION/'

nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

domsize=768

nc_data3D = Dataset(nc_in3D)
nc_data2D = Dataset(nc_in2D)
varis3D = nc_data3D.variables
varis2D = nc_data2D.variables

print 'loading 2D variables'
nc_data2D = Dataset(nc_in2D)
varis2D = nc_data2D.variables

PW = varis2D['PW'][:]

print 'loading 3D variables'
nc_data3D = Dataset(nc_in3D)
varis3D = nc_data3D.variables

t3D = varis3D['time'][:]
t2D = varis2D['time'][:]
x = varis2D['x'][:]
y = varis2D['y'][:]
z = varis3D['z'][:]
p = varis3D['p'][:]
p = p*1e2

#averaging time period for 2D & 3D fields
ntave2D=24
ntave3D=4

nt3D = t3D.size
nt2D = t2D.size
nz = z.size
nx = x.size
ny = y.size

xx, yy = np.meshgrid(x, y)
times = np.arange(t3D[0], np.max(t3D))

#2D fields
ntrunc = PW.shape[0]%ntave2D
PW = PW[ntrunc:,:,:]
#PW_tmp = PW.reshape(nt2D/ntave2D, ntave2D, nx, ny)
#PW_tave = np.mean(PW_tmp, axis=1)

#EDIT: time period to look at (t in days, nave averaging period in days)
t=-35
nave=5

nave2D=5*ntave2D
nave3D=5*ntave3D
t2=t*ntave2D
t3=t*ntave3D
#number of blocks to find PW max
db=1
print 'calulating max PW, blocked'
#PW_t = np.mean(PW_tave[t-nave:t,:,:], axis=0)
PW_tave = np.mean(PW[t2-nave2D:t2,:,:], axis=0)
PW_blocked = blockave2D(PW_tave, db)
PWxy_sorted = blockxysort2D(PW_tave, xx, yy, db)

PWsort = PWxy_sorted.keys()
PWxycoords = PWxy_sorted.values()

mcenter = PWxycoords[-1]

nbins=[100, z[:-1]]

buoyflux = np.load(fin + 'buoyflux_evolution.npy')

print 'calculating temporal averages'
#average over nave days 
buoyflux_tave = np.mean(buoyflux[t-nave:t,:,:,:], axis=0)

#fieldmeans has shape (redges.size, zedges.size) = nbins
print '2d contouring'
redges, zedges, fieldmeans = radprof3D(buoyflux_tave, xx, yy, z[:-1], mcenter, nbins=nbins)
rbin_centers = (redges[1:] + redges[:-1])/2.
zbin_centers = zedges

#rr has shape (zedges.size, redges.size) = nbins.T
rr, zz = np.meshgrid(rbin_centers, zbin_centers)

plt.figure()
ax=plt.gcf().gca()
vmin=np.ma.masked_invalid(fieldmeans).min()
vmax=np.ma.masked_invalid(fieldmeans).max()
#EDIT: COLORBAR FORMATTING. IF NOT MONOTONIC, SWITCH TO DIVERGING COLORBAR, OTHERWISE
#USE SEQUENTIAL COLORBAR. SET COLOR MAPPING LIMITS MANUALLY. 
if (vmin < 0 and vmax > 0):
    if np.abs(vmin) < np.abs(vmax):
        vmin = vmin - (vmax - np.abs(vmin))
    else:
        vmax = vmax + (np.abs(vmin) - vmax)
    cmap = cm.RdBu_r
else:
    cmap = cm.YlGnBu
print 'plotting'
fieldmeans = np.transpose(fieldmeans)
plt.pcolormesh(rr/(1e3*domsize), zz/(1e3), fieldmeans, vmin=vmin, vmax=vmax, cmap=cmap)
plt.xlabel('fractional distance from moist region center, relative to domain size')
plt.ylabel('z (km)')
cb=plt.colorbar()
cb.set_label(r'W/m$^3$')
ax.set_ylim(0,26)
plt.title('buoyancy flux, (W/m$^3$), day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$, x bin width = {:2.2f}'.format(t3D[t3-nave3D], t3D[t3], domsize, 1./nbins[0]))
plt.savefig(fout + 'buoyfluxradialxsection_day{:2.0f}to{:2.0f}.pdf'.format(t3D[t3-nave3D], t3D[t3]))
plt.close()




