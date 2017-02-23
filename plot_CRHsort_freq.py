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
fout = '/Users/cpatrizio/figures/SST302/SAM_aggr130days_768km_64vert_ubarzero_CRHSORT/'

nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

domsize=768

print 'loading variables'
#t3D = varis3D['time'][:]

nc_data2D = Dataset(nc_in2D)
varis2D = nc_data2D.variables
t2D = varis2D['time'][:]
x = varis2D['x'][:]
y = varis2D['y'][:]
#z = varis3D['z'][:]
#p = varis3D['p'][:]
#p = p*1e2

#nt3D = t3D.size
nt2D = t2D.size
#nz = z.size
nx = x.size
ny = y.size

#averaging time period for 2D & 3D fields
ntave2D=24
ntave3D=4

print 'loading netCDF files'

nc_fs = glob.glob(fpath2D + '*256x256*3000m*302K.nc')

CRHs = np.array([])

for nc_in2D in nc_fs:
#nc_data3D = Dataset(nc_in3D)
    nc_data2D = Dataset(nc_in2D)
    #varis3D = nc_data3D.variables
    varis2D = nc_data2D.variables
    PW = varis2D['PW'][:]
    SWVP = varis2D['SWVP'][:]
    CRH = PW/SWVP
    CRH_tmp = CRH.reshape(nt2D/ntave2D, ntave2D, nx, ny)
    CRH_tave = np.mean(CRH_tmp, axis=1)
    db=16
    CRHs_temp = blockave3D(CRH_tave, db)
    #        buoyflux = np.vstack((buoyflux, buoyflux_temp)) if buoyflux.size else buoyflux_temp
    CRHs = np.vstack((CRHs, CRHs_temp)) if CRHs.size else CRHs_temp
    ntimes = CRHs.shape[0]



    

#calculate CRH
PW = varis2D['PW'][:]
SWVP = varis2D['SWVP'][:]
CRH = PW/SWVP
CRH_tmp = CRH.reshape(nt2D/ntave2D, ntave2D, nx, ny)

#daily average
CRH_tave = np.mean(CRH_tmp, axis=1)

db=16
nblocks = (nx/db)*(ny/db)
CRHranks = np.arange(nblocks)

#calculate 2D (CRH rank, time) map of CRH frequency  ###
times = np.arange(0, 130)
CRHs =  np.zeros((times.size, nblocks))

#do block averaging
CRHs = blockave3D(CRH_tave, db)
#for i, t in enumerate(times):
#   CRHs[i,:] = np.sort(blockave2D(CRH_tave[i,:,:], db).flatten())
    
CRHrankss, tt = np.meshgrid(CRHranks, times)
tcoords = tt.flatten()
CRHcoords = CRHs.flatten()
CRHfreqs, xedges, yedges = np.histogram2d(CRHcoords, np.transpose(tcoords), bins=(256,130))

extent = [CRHranks[0], CRHranks[-1], yedges[-1], yedges[0]]

#plot temporal evolution of CRH frequency distribution
plt.figure(3)
plt.imshow(np.transpose(CRHfreqs), extent=extent, vmin=0, vmax=40, cmap=cm.bone_r)
cb = plt.colorbar()
cb.set_label('frequency')
plt.title(r'temporal evolution of CRH frequency distribution, domain size {:3.0f} km$^2$'.format(domsize))
plt.xlabel('CRH rank')
plt.ylabel('time (days)')
plt.savefig(fout + 'CRHsort_freqvst.pdf')
plt.close()
t1=5
t2=110
#plot frequency distribution of CRH at specific times (t1, t2)
plt.figure(4)
plt.hist(CRHs[t1,:], bins=np.linspace(0,1,30), histtype='stepfilled', label='day {:2.1f}'.format(t1))
plt.hist(CRHs[t2,:], bins=np.linspace(0,1,30), histtype='stepfilled', label='day {:2.1f}'.format(t2))
plt.xlabel('CRH')
plt.ylabel('frequency')
plt.legend()
plt.savefig(fout + 'CRHsort_freqday{:2.1f}and{:2.1f}.pdf'.format(t1, t2))
plt.close()

