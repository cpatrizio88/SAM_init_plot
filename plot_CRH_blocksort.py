from netCDF4 import Dataset
import site
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib
from thermolib.wsat import wsat
from SAM_init_plot.block_fns import blockave2D

rho_w = 1000 #density of water
g=9.81 #gravitational acceleration
p_s = 1015

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/figures/SST302/SAM_aggrday90to130_768km_64vert_ubarzero_3D/'

nc_in3D = glob.glob(fpath3D + '*256x256*3000m*day90to130*302K.nc')[0]
nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

T_s = 302

nc_data3D = Dataset(nc_in3D)
nc_data2D = Dataset(nc_in2D)
varis3D = nc_data3D.variables
varis2D = nc_data2D.variables

ti=-1
ntave=4

x = varis3D['x'][:]
y = varis3D['y'][:]
xx, yy = np.meshgrid(x, y)
z = varis3D['z'][:]
p = varis3D['p'][:]
p = p*1e2
qv = varis3D['QV'][:]
qv = qv*1e-3
T = varis3D['TABS'][:]
#T_tave = np.mean(T[ti-ntave,:,:,:])
#qv_tave = np.mean(qv[ti-ntave,:,:,:])
delp3D = np.ones((x.size, y.size, z.size-1))
delp = -np.diff(p)
delp3D[:,:,:] = delp
delp3D = np.transpose(delp3D)

PW = varis2D['PW'][:]
PW_tave = np.mean(PW[ti-ntave:ti,:,:], axis=0)

#PW = 1/(g*rho_w)*np.sum(np.multiply(delp3D, qv_t[:-1,:,:]), axis=0)
#SPW = np.zeros(PW.shape)
#
#for i, plev in enumerate(p[:-1]):
#    print i
#    SPW[:,:] = 1/(g*rho_w)*np.sum(np.multiply(delp3D, wsat(T_i[i,:,:], plev)))
#CRH = PW/SPW

#width of block in units of dx
db=16
PW_blocked = blockave2D(PW_tave, db)

#plt.figure(1)
#plt.contour(xx/1e3, yy/1e3, PW_tave, 20, colors='k', alpha=0.5)
#plt.contourf(xx/1e3, yy/1e3, PW_tave, 20, cmap=cm.RdYlBu_r, zorder=0)
#plt.title('non-blocked')
#plt.show()
#
#plt.figure(2)
#plt.contour(xx[::db,::db]/1e3, yy[::db,::db]/1e3, PW_blocked, 20, colors='k', alpha=0.5)
#plt.contourf(xx[::db,::db]/1e3, yy[::db,::db]/1e3, PW_blocked, 20, cmap=cm.RdYlBu_r, zorder=0)
#plt.title('blocked')
#plt.show()





