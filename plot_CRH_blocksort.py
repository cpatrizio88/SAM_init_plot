from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker
import matplotlib
from CRH_blocksort import CRH_blocksort
from wsat import wsat



rho_w = 1000 #density of water
g=9.81 #gravitational acceleration
p_s = 1015

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fout = '/Users/cpatrizio/figures/SST302/SAM_aggrday90to130_768km_64vert_ubarzero_3D/'

nc_in = glob.glob(fpath + '*3000m*day90to130*302K.nc')[0]

nc_data= Dataset(nc_in)
varis = nc_data.variables

delt=6
ti=-1

x = varis['x'][:]
y = varis['y'][:]
z = varis['z'][:]
p = varis['p'][:]
p = np.concatenate(([p_s], p))
p = p*1e2
qv = varis['QV'][:]
T = varis['TABS'][:]
T_t = T[ti,:,:,:]
qv_t = qv[ti,:,:,:]
delp3D = np.ones((x.size, y.size, z.size))
delp = -np.diff(p)
delp3D[:,:,:] = delp
delp3D = np.transpose(delp3D)
    
PW = 1/(g*rho_w)*np.sum(np.multiply(delp3D, qv_t), axis=0)

#CRH = CRH_blocksort(varis, 'null', 0, -1, 0)

