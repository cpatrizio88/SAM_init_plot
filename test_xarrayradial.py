import xarray as xr
import numpy as np
import glob
from SAM_init_plot.block_fns import blockave2D, blockave3D, blockxysort2D


fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'

nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

nave3D=4
nave2D=24

ds3D = xr.open_dataset(nc_in3D, chunks={'time': nave3D})
ds2D = xr.open_dataset(nc_in2D, chunks={'time': nave2D})

t3D = ds3D['time'][:]
t2D = ds2D['time'][:]
x = ds3D['x'][:]
y = ds3D['y'][:]
z = ds3D['z'][:]
p = ds3D['p'][:]
p = p*1e2

#varnames = ['QV', 'W', 'U', 'QRAD', 'QN', 'QP']
#varnames=['QV']
varnames=['QRAD']

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

PW = ds2D['PW'][:]

#2D fields
ntrunc = PW.shape[0]%ntave2D
PW = PW[ntrunc:,:,:]
#PW_tmp = PW.reshape(nt2D/ntave2D, ntave2D, nx, ny)
#PW_tave = np.mean(PW_tmp, axis=1)

#time period to look at
t=-25*ntave3D
nave=5
nave2D=5*ntave2D
nave3D=5*ntave3D
#number of blocks to find PW max
db=1
print 'calulating max PW, blocked'
#PW_t = np.mean(PW_tave[t-nave:t,:,:], axis=0)
PW_t = np.mean(PW[t-nave2D:t,:,:], axis=0)
PW_blocked = blockave2D(PW_t, db)
PWxy_sorted = blockxysort2D(PW_t, xx, yy, db)

PWsort = PWxy_sorted.keys()
PWxycoords = PWxy_sorted.values()

mcenter = PWxycoords[-1]

nbins=[100, z]
