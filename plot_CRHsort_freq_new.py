from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from thermolib.wsat import wsat
from thermolib.constants import constants
import SAM_init_plot.block_fns
from SAM_init_plot.block_fns import blockave2D, blockave3D

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggr150days_64vert_ubarzero_CRHSORT/'

#nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day90to130*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*512x512*3000m*day90to140*302K.nc')[0]


#EDIT 
domsize=3072

#EDIT 
nc_fs = glob.glob(fpath2D + '*1024x1024*3000m*302K.nc')

#EDIT 
tinit = 90

print 'loading variables'
#t3D = varis3D['time'][:]

nc_data2D = Dataset(nc_fs[0])
varis2D = nc_data2D.variables
x = varis2D['x'][:]
y = varis2D['y'][:]
#z = varis3D['z'][:]
#p = varis3D['p'][:]
#p = p*1e2

#nt3D = t3D.size
#nt2D = t2D.size
#nz = z.size
nx = x.size
ny = y.size

#averaging time period for 2D & 3D fields
ntave2D=24
ntave3D=4

print 'loading netCDF files'

CRHs = np.array([])
PWs = np.array([])

for nc_in2D in nc_fs:
#nc_data3D = Dataset(nc_in3D)
    print 'loading', nc_in2D
    nc_data2D = Dataset(nc_in2D)
    #varis3D = nc_data3D.variables
    varis2D = nc_data2D.variables
    t2D = varis2D['time'][:]
    nt2D = t2D.size
    trunc = nt2D%ntave2D
    PW = varis2D['PW'][trunc:]
    #SWVP = varis2D['SWVP'][trunc:]
    #CRH = PW/SWVP
    #CRH_tmp = CRH.reshape(nt2D/ntave2D, ntave2D, nx, ny)
    PW_tmp = PW.reshape(nt2D/ntave2D, ntave2D, nx, ny)
    PW_tave = np.mean(PW_tmp, axis=1)
    #CRH_tave = np.mean(CRH_tmp, axis=1)
    db=16
    #CRHs_temp = blockave3D(CRH_tave, db)
    PWs_temp = blockave3D(PW_tave, db)
    #        buoyflux = np.vstack((buoyflux, buoyflux_temp)) if buoyflux.size else buoyflux_temp
    #CRHs = np.vstack((CRHs, CRHs_temp)) if CRHs.size else CRHs_temp
    PWs = np.vstack((PWs, PWs_temp)) if PWs.size else PWs_temp

#ntimes = CRHs.shape[0]
ntimes = PWs.shape[0]

times = np.arange(tinit, tinit+ntimes)
nblocks = (nx/db)*(ny/db)
#CRHranks = np.arange(nblocks)

#CRHbins = np.linspace(0, 1, nblocks)
PWbins = np.linspace(0, 85, nblocks)

#CRHrankss, tt = np.meshgrid(CRHbins, times)
PWbinss, tt = np.meshgrid(PWbins, times)
tcoords = tt.flatten()
#CRHcoords = CRHs.flatten()
PWcoords = PWs.flatten()
#CRHfreqs, xedges, yedges = np.histogram2d(CRHcoords, np.transpose(tcoords), bins=(CRHbins, ntimes))
PWfreqs, xedges, yedges = np.histogram2d(PWcoords, np.transpose(tcoords), bins=(PWbins, ntimes))

extent = [xedges[0], xedges[-1], yedges[-1], yedges[0]]

#plot temporal evolution of CRH frequency distribution
plt.figure()
#plt.imshow(np.transpose(CRHfreqs), extent=extent, vmin=0, vmax=15, cmap=cm.magma_r, aspect=0.5)
plt.imshow(np.transpose(PWfreqs), extent=extent, vmin=0, vmax=10, cmap=cm.magma_r, aspect=1)
cb = plt.colorbar()
cb.set_label('frequency')
plt.title(r'temporal evolution of PW frequency distribution, domain size ({:3.0f} km)$^2$, block size ({:3.0f} km)$^2$'.format(domsize, (db*np.diff(x)[0])/1e3))
plt.xlabel('PW (mm)')
plt.ylabel('time (days)')
plt.savefig(fout + 'PWsort_day{:3.0f}to{:3.0f}_freq.pdf'.format(times[0], times[-1]))
plt.close()
t1=5
t2=110
#plot frequency distribution of CRH at specific times (t1, t2)
#plt.figure(4)
#plt.hist(CRHs[t1,:], bins=np.linspace(0,1,30), histtype='stepfilled', label='day {:2.1f}'.format(t1))
#plt.hist(CRHs[t2,:], bins=np.linspace(0,1,30), histtype='stepfilled', label='day {:2.1f}'.format(t2))
#plt.xlabel('CRH')
#plt.ylabel('frequency')
#plt.legend()
#plt.savefig(fout + 'CRHsort_freqday{:2.1f}and{:2.1f}.pdf'.format(t1, t2))
#plt.close()

