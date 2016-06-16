from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from thermolib.wsat import wsat
from thermolib.constants import constants
import SAM_init_plot.block_fns
from SAM_init_plot.block_fns import blockave1D

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/figures/SST302/SAM_aggr100days_12288km_64vert_ubarzero_CRHSORT/'

nc_in3D = glob.glob(fpath3D + '*4096*64*3000m*100days*302K.nc')[0]
nc_in2D = glob.glob(fpath2D + '*4096*64*3000m*100days*302K.nc')[0]

domsize=12288

print 'loading netCDF files'
nc_data3D = Dataset(nc_in3D)
nc_data2D = Dataset(nc_in2D)
varis3D = nc_data3D.variables
varis2D = nc_data2D.variables

print 'loading variables'
t3D = varis3D['time'][:]
t2D = varis2D['time'][:]
x = varis3D['x'][:]
z = varis3D['z'][:]
p = varis3D['p'][:]
p = p*1e2

nt3D = t3D.size
nt2D = t2D.size
nz = z.size
nx = x.size

#averaging time period for 2D & 3D fields
ntave2D=24
ntave3D=8

#calculate CRH
PW = varis2D['PW'][:]
SWVP = varis2D['SWVP'][:]
CRH = PW/SWVP
CRH_tmp = CRH.reshape(nt2D/ntave2D, ntave2D, nx)

#daily average
CRH_tave = np.mean(CRH_tmp, axis=1)

db=16
nblocks = (nx/db)
CRHranks = np.arange(nblocks)

#calculate 2D (CRH rank, time) map of CRH frequency  ###
times = np.arange(0, t3D[-1])
CRHs =  np.zeros((times.size, nblocks))

#do block averaging
for i, t in enumerate(times):
   CRHs[i,:] = np.sort(blockave1D(CRH_tave[i,:], db).flatten())
    
CRHrankss, tt = np.meshgrid(CRHranks, times)
tcoords = tt.flatten()
CRHcoords = CRHs.flatten()
CRHfreqs, xedges, yedges = np.histogram2d(CRHcoords, np.transpose(tcoords), bins=(256,100))

extent = [CRHranks[0], CRHranks[-1], yedges[-1], yedges[0]]

#plot temporal evolution of CRH frequency distribution
plt.figure(1)
plt.imshow(np.transpose(CRHfreqs), extent=extent, vmin=0, vmax=15, cmap=cm.bone_r)
cb = plt.colorbar()
cb.set_label('frequency')
plt.title('temporal evolution of CRH frequency distribution, domain size {:5.0f} km'.format(domsize))
plt.xlabel('CRH rank')
plt.ylabel('time (days)')
plt.savefig(fout + 'CRHsort_freqvst.pdf')
plt.close()

t1=10
t2=-10
#plot frequency distribution of CRH at specific times (t1, t2)
plt.figure(2)
plt.hist(CRHs[t1,:], bins=np.linspace(0,1,30), histtype='stepfilled', label='day {:2.1f}'.format(times[t1]))
plt.hist(CRHs[t2,:], bins=np.linspace(0,1,30), histtype='stepfilled', label='day {:2.1f}'.format(times[t2]))
plt.xlabel('CRH')
plt.ylabel('frequency')
plt.title('frequency distribution of CRH, domain size {:5.0f} km'.format(domsize))
plt.legend()
plt.savefig(fout + 'CRHsort_freqday{:2.1f}and{:2.1f}.pdf'.format(times[t1], times[t2]))
plt.close()

