from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
import SAM_init_plot.block_fns
import SAM_init_plot.misc_fns
from SAM_init_plot.misc_fns import raddist, radprof
from SAM_init_plot.block_fns import blockave2D, blockxysort2D, xysort

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/figures/SST302/SAM_aggrday130_768km_64vert_ubarzero_RADIAL/'

nc_in2D = glob.glob(fpath + '*256x256*3000m*130days*302K.nc')[0]

print 'loading 2D variables'
nc_data2D = Dataset(nc_in2D)
varis2D = nc_data2D.variables

PW = varis2D['PW'][:]
t2D = varis2D['time'][:]
x = varis2D['x'][:]
y = varis2D['y'][:]
nt2D = t2D.size
nx = x.size
ny = y.size

ntave2D=24

PW_tmp = PW.reshape(nt2D/ntave2D, ntave2D, nx, ny)
PW_tave = np.mean(PW_tmp, axis=1)

xx, yy = np.meshgrid(x, y)

t=-20
db=16

PW_t = PW_tave[t,:,:]

PW_blocked = blockave2D(PW_t, db)

PWxy_sorted = blockxysort2D(PW_t, xx, yy, db)

PWsort = PWxy_sorted.keys()
xycoords = PWxy_sorted.values()

mcenter = xycoords[-1]

d = raddist(xx, yy, mcenter)

#put sum of PW in bins of r
#nr is r bins (indices give the value of r)
#PWrsums is sum of PW in bins of r
nr, PWrsums = radprof(PW_t, xx, yy, mcenter)

PWrprof = PWrsums/nr
rs = np.arange(0, len(nr))

plt.figure(1)
plt.contourf(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, PW_blocked, 30, cmap=cm.RdYlBu_r)
plt.colorbar()
plt.title('PW (mm)')
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.show()

plt.figure(2)
plt.contourf(xx/1e3, yy/1e3, d/1e3, 30, cmap=cm.RdYlBu_r)
plt.title('distances')
plt.show()

plt.figure(3)
plt.plot(rs/1e3, PWrprof, 'k,')
plt.xlabel('distance from moist region center (km)')
plt.ylabel('PW (mm)')
plt.show()



