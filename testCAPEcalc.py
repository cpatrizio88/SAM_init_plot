from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker
import matplotlib
from SAM_init_plot.block_fns import blockave3D, blockave2D
from thermolib.constants import constants
from thermolib.constants import constants
import thermolib.thermo as thermo
from thermolib.findTmoist import findTmoist
from thermolib.wsat import wsat

c = constants()

matplotlib.rcParams.update({'font.size': 26})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 2})
matplotlib.rcParams.update({'legend.fontsize': 22})

plt.style.use('seaborn-white')
fpathSTAT = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
fpath2D =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fpath3D = '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
foutdata = '/Users/cpatrizio/data/SST302/'
#fout =  '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_MAPSNEW/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to140_64vert_ubarzero_MAPSNEW/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to150_64vert_ubarzero_MAPSNEW/'


#nc_STAT = glob.glob(fpathSTAT + '*256x256*3000m*250days*302K.nc')[0]
#nc_in = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
 
nc_STAT = glob.glob(fpathSTAT + '*512x512*3000m*195days*302K.nc')[0]
nc_in = glob.glob(fpath2D + '*512x512*3000m*day180to195*302K.nc')[0]
nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day180to195*302K.nc')[0]

#nc_STAT = glob.glob(fpathSTAT + '*1024x1024*3000m*170days*302K.nc')[0]
#nc_in = glob.glob(fpath2D + '*1024x1024*3000m*day090to130*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]

#domsize=768
domsize=1536
#domsize=3072

nc_data= Dataset(nc_in)
nc_dataSTAT = Dataset(nc_STAT)
nc_data3D = Dataset(nc_in3D)
varis = nc_data.variables
varis3D = nc_data3D.variables
varisSTAT = nc_dataSTAT.variables
times2D = varis['time'][:]

#one moist adiabat, and subtract from temperature field

p = varis3D['p'][:]
z = varis3D['z'][:]
x = varis['x'][:]
y = varis['y'][:]

xx, yy = np.meshgrid(x, y)

nz = z.size
nx = x.size
ny = y.size


t = -1

p_s = p[0]
T_s = varisSTAT['TABS'][t,0]
q_sat = wsat(T_s, p_s*1e2)

Tenv = varis3D['TABS'][t,:,:,:]

QV = varis3D['QV'][t,:,:,:]
QV = QV*1e-3

Tenv_ave = np.mean(np.mean(Tenv, axis=2), axis=1)


thetae0 = thermo.theta_e(T_s, p_s*1e2, q_sat, 0) #theta_e in moist region 

Tadiabat = findTmoist(thetae0, p*1e2)

qvadiabat = []

for i, plev in enumerate(p*1e2):
    Tadb = Tadiabat[i]
    qvadiabat = qvadiabat + [wsat(Tadb, plev)]
    
qvadiabat = np.array(qvadiabat)
    
TVenv = Tenv*(1+0.61*QV)

TVenv_ave = np.mean(np.mean(TVenv, axis=2), axis=1)
    
TVadiabat = Tadiabat*(1+0.61*qvadiabat)


TVadiabat3D = np.zeros((ny, nx, nz))
TVadiabat3D[:,:,:] = TVadiabat
TVadiabat3D = TVadiabat3D.T

Tpert = TVadiabat3D - TVenv

Tpert_ave = np.mean(np.mean(Tpert, axis=2), axis=1)

znb_i = np.where(Tpert_ave[::-1] > 0)[0][0]
z_nb = z[::-1][znb_i]
znb_i = np.where(z > z_nb)[0][0]
zL_i = np.where(Tpert_ave > 0)[0][0]

zL_i = 0

delz = np.diff(z)
delz3D = np.zeros((ny, nx, nz-1))
delz3D[:,:,:] = delz
delz3D = delz3D.T
        
CAPE = c.g*np.sum(np.multiply(Tpert[zL_i:znb_i,:,:]/TVenv[zL_i:znb_i,:,:], delz3D[zL_i:znb_i,:,:]), axis=0)

#CIN = c.g*np.sum(np.multiply(Tpert[:zL_i,:,:]/TVenv[:zL_i,:,:], delz3D[:zL_i,:,:]), axis=0)
#zf_i = zf_i[zf_i > 20




plt.figure(1)
plt.plot(Tadiabat[:znb_i], z[:znb_i]/1e3, 'k', label = r'moist adiabat, $\theta_e$ = {:3.1f} K'.format(thetae0))
plt.plot(Tenv_ave[:znb_i], z[:znb_i]/1e3, 'b', label = 'domain-mean temperature')
plt.xlabel('temperature (K)')
plt.title('temperature profile on day {:3.0f}, domain width = {:4.0f} km'.format(times2D[t], domsize))
plt.ylabel('z (km)')
plt.legend()
plt.show()


plt.figure(2)
plt.plot(TVadiabat, z/1e3, 'k', label = r'moist adiabat, $\theta_e$ = {:3.1f} K'.format(thetae0))
plt.plot(TVenv_ave, z/1e3, 'b', label = 'domain-mean virtual temperature')
plt.xlabel('virtual temperature (K)')
plt.title('virtual temperature profile on day {:3.0f}, domain width = {:4.0f} km'.format(times2D[t], domsize))
plt.ylabel('z (km)')
plt.legend()
plt.show()

vals = np.arange(300, 1200, 20)

plt.figure(3)
plt.contourf(xx/1e3, yy/1e3, CAPE, vals, cmap=cm.RdYlBu_r, zorder=0)
cb = plt.colorbar()
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.title('CAPE')
plt.show()

#vals = np.arange(-10, 0, 0.1)
#
#plt.figure(4)
#plt.contourf(xx/1e3, yy/1e3, CIN, vals, cmap=cm.RdYlBu, zorder=0)
#cb = plt.colorbar()
#plt.xlabel('x (km)')
#plt.ylabel('y (km)')
#plt.title('CIN')
#plt.show()










