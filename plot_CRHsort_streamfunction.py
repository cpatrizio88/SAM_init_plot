from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from thermolib.wsat import wsat
from thermolib.constants import constants
import SAM_init_plot.block_fns
from SAM_init_plot.block_fns import blocksort3D

c = constants()

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/figures/SST302/SAM_aggrday90to130_768km_64vert_ubarzero_CRHSORT/'

nc_in3D = glob.glob(fpath3D + '*256x256*3000m*day90to130*302K.nc')[0]
nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

domsize=768

print 'loading netCDF files'
nc_data3D = Dataset(nc_in3D)
nc_data2D = Dataset(nc_in2D)
varis3D = nc_data3D.variables
varis2D = nc_data2D.variables

print 'loading variables'
t3D = varis3D['time'][:]
t2D = varis2D['time'][:]
x = varis3D['x'][:]
y = varis3D['y'][:]
z = varis3D['z'][:]
p = varis3D['p'][:]
p = p*1e2

nt3D = t3D.size
nt2D = t2D.size
nz = z.size
nx = x.size
ny = y.size

#calculate CRH
PW = varis2D['PW'][:]
SWVP = varis2D['SWVP'][:]
CRH = PW/SWVP

#need w to calculate streamfuncion
w = varis3D['W'][:]

#need T to calculate density
T = varis3D['TABS'][:]

rhoz = p/(c.Rd*np.mean(np.mean(T, axis=3), axis=2))

#averaging time period for 2D & 3D fields
ntave2D=24
ntave3D=4

print 'calculating daily averages'
#reshape to do daily averages
T_tmp = T.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
w_tmp = w.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
rho_tmp = rhoz.reshape(nt3D/ntave3D, ntave3D, nz)
CRH_tmp = CRH.reshape(nt2D/ntave2D, ntave2D, nx, ny)

#daily average
CRH_tave = np.mean(CRH_tmp, axis=1)
T_tave = np.mean(T_tmp, axis=1)
w_tave = np.mean(w_tmp, axis=1)
rhoz_tave = np.mean(rho_tmp, axis=1)

#times in days
times = np.arange(int(t3D[0]), int(t3D[-1]))


#number of grid points in a block
db=16

###### USER EDIT: time to look at streamfunction
t = -30

#sort CRH and w
CRH_wsort = blocksort3D(CRH_tave[t,:,:], w_tave[t,:,:,:], db)
CRHsort = CRH_wsort[0]
wsort = CRH_wsort[1]

#calculate the streamfunction
print 'calculating streamfunction'
psi = np.zeros((z.size, CRHsort.size))
for i in range(z.size):
    for j in range(1, CRHsort.size):
      psi[i,j] = psi[i, j-1] + np.sum(rhoz_tave[t, i]*wsort[i,:j])
      
CRHranks = np.arange(np.size(CRHsort))
CRHrankss, zz = np.meshgrid(CRHranks, z)

########USER EDIT: variable to overlay with streamfunction
varname = 'QV'
vari = varis3D[varname]
var_tmp = vari[:].reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
var_tave = np.mean(var_tmp, axis=1)

CRH_varsort = blocksort3D(CRH_tave[t,:,:], var_tave[t,:,:,:], db)
varsort = CRH_varsort[1]

plt.figure(1)
ax=plt.gcf().gca()
plt.contour(CRHrankss, zz/1e3, psi, 15)
plt.contourf(CRHrankss, zz/1e3, varsort, 20, cmap=cm.RdYlBu_r)
ax.set_ylim((0,15))
plt.xlabel('CRH rank')
plt.ylabel('height (km)')
plt.title(r'streamfunction at day {:2.1f}, domain size {:3.0f} km$^2$'.format(times[t], domsize))
cb = plt.colorbar()
cb.set_label('{:s} ({:3s})'.format(varname, vari.units.strip()))
plt.savefig(fout + 'CRHsort_streamfnday{:2.1f}_{:s}.pdf'.format(times[t], varname))
plt.close()






      
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    