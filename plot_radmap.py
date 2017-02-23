from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker
import matplotlib
from thermolib.constants import constants
from SAM_init_plot.misc_fns import raddist, radprof
from SAM_init_plot.block_fns import blockave2D, blockxysort2D, xysort

c = constants()

fpathSTAT = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
fpath2D =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fpath3D = '/Users/cpatrizio/SAM6.10.8/OUT_3D/'

fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr250days__64vert_ubarzero_RADDISTMAPS/'

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 2})

nc_STAT = glob.glob(fpathSTAT + '*256x256*3000m*250days*302K.nc')[0]
nc_in = glob.glob(fpath2D + '*256x256*3000m*day230to250*302K.nc')[0]
nc_in3D = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]

#nc_STAT = glob.glob(fpathSTAT + '*512x512*3000m*170days*302K.nc')[0]
#nc_in = glob.glob(fpath2D + '*512x512*3000m*day140to170*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day155to170*302K.nc')[0]

#nc_STAT = glob.glob(fpathSTAT + '*1024x1024*3000m*150days*302K.nc')[0]
#nc_in = glob.glob(fpath2D + '*1024x1024*3000m*day140to150*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*1024x1024*3000m*day140to150*302K.nc')[0]

domsize=768
#domsize=1536
#domsize=3072

nc_data= Dataset(nc_in)
nc_dataSTAT = Dataset(nc_STAT)
nc_data3D = Dataset(nc_in3D)
varis2D = nc_data.variables
varis3D = nc_data3D.variables
varisSTAT = nc_dataSTAT.variables
times2D = varis2D['time'][:]



x = varis2D['x'][:]
y = varis2D['y'][:]


db=1


xx, yy = np.meshgrid(x, y)

nave=5
ntave3D=4
ntave2D=24
t=-1

aveperiod2D = nave*ntave2D
aveperiod3D = nave*ntave3D

PW = varis2D['PW'][t-aveperiod2D:t,:,:]

W = varis3D['W'][t-aveperiod3D:t,:,:,:]

W_tave = np.mean(W, axis=0)
W_ztave = np.mean(W_tave, axis=0)


PW_tave = np.mean(PW, axis=0)

db=1
dbw=2

W_ztave = blockave2D(W_ztave, dbw)


PW_blocked = blockave2D(PW_tave, db)

PWxy_sorted = blockxysort2D(PW_tave, xx, yy, db)
#fieldxy_sorted = blockxysort2D(field_t, xx, yy, db)

PWsort = PWxy_sorted.keys()
PWxycoords = PWxy_sorted.values()
#fieldsort = fieldxy_sorted.keys()
#fieldxycoords = fieldxy_sorted.values()

mcenter = PWxycoords[-1]

d = raddist(xx, yy, mcenter)

rvals = np.arange(0, 1./np.sqrt(2), 5./domsize)
PWvals = np.arange(0, 90, 5)

Wvals = np.arange(-0.004, 0, 0.0001)


plt.figure(1)
plt.contourf(xx/(domsize*1e3), yy/(domsize*1e3), d/(domsize*1e3), levels=rvals, cmap=cm.YlGnBu, zorder=0)
plt.plot(mcenter[0]/(domsize*1e3), mcenter[1]/(domsize*1e3), 'x', mew=3, zorder=1)
plt.xlabel('x/L')
plt.ylabel('y/L')
plt.title('radial distance from moist region center (km), domain size (L) = {:d} km'.format(domsize))
plt.colorbar()
plt.savefig(fout + '{:d}km_rdistmap_day{:3.0f}to{:3.0f}'.format(domsize, times2D[t-aveperiod2D], times2D[t]))
plt.close()

plt.figure(2)
plt.contourf(xx/(domsize*1e3), yy/(domsize*1e3), PW_tave, levels=PWvals, cmap=cm.YlGnBu, zorder=0)
plt.plot(mcenter[0]/(domsize*1e3), mcenter[1]/(domsize*1e3), 'x', mew=3, zorder=1)
plt.xlabel('x/L')
plt.ylabel('y/L')
plt.title('PW (mm), domain size (L) = {:d} km'.format(domsize))
plt.colorbar()
plt.savefig(fout + '{:d}km_PWmap_day{:3.0f}to{:3.0f}'.format(domsize, times2D[t-aveperiod2D], times2D[t]))
plt.close()

plt.figure(3)
plt.contourf(xx[::dbw,::dbw]/(domsize*1e3), yy[::dbw, ::dbw]/(domsize*1e3), W_ztave, levels=Wvals, cmap=cm.YlGnBu, zorder=0)
plt.plot(mcenter[0]/(domsize*1e3), mcenter[1]/(domsize*1e3), 'x', mew=3, zorder=1)
plt.xlabel('x/L')
plt.ylabel('y/L')
plt.title(r'vertically-averaged W (m/s), block-averaged over ({:2.0f} km)$^2$ domain size (L) = {:d} km'.format((dbw*np.diff(x)[0])/1e3, domsize))
plt.colorbar()
plt.savefig(fout + '{:d}km_vertaveWmap_day{:3.0f}to{:3.0f}'.format(domsize, times2D[t-aveperiod2D], times2D[t]))
plt.close()

