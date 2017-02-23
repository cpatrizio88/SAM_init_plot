from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from thermolib.wsat import wsat
from thermolib.constants import constants
import SAM_init_plot.block_fns
from SAM_init_plot.block_fns import blockave2D, blocksort2D, blocksort3D, vertint

c = constants()

p_s = 1015
T_s = 302

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (18, 12)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_CRHSORT/'

nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
nc_in2D = glob.glob(fpath2D + '*256*256*3000m*130days*302K.nc')[0]

domsize=768

nc_data3D = Dataset(nc_in3D)
nc_data2D = Dataset(nc_in2D)
varis3D = nc_data3D.variables
varis2D = nc_data2D.variables

t3D = varis3D['time'][:]
t2D = varis2D['time'][:]
x = varis3D['x'][:]
y = varis3D['y'][:]
z = varis3D['z'][:]
p = varis3D['p'][:]
p = p*1e2

#2D variables
PW = varis2D['PW'][:]
SWVP = varis2D['SWVP'][:]

CRH = PW/SWVP

#3D variables
print 'loading 3D variables'
qv = varis3D['QV'][:]
qv = qv*1e-3
T = varis3D['TABS'][:]
w = varis3D['W'][:]
Qr = varis3D['QRAD'][:]

#1D variables
rhoz = p/(c.Rd*np.mean(np.mean(T, axis=3), axis=2))

#averaging time period for 2D & 3D fields
ntave2D=24
ntave3D=4

nt3D = t3D.size
nt2D = t2D.size
nz = z.size
nx = x.size
ny = y.size

dt3D = np.diff(t3D)[0]

dt3D = dt3D*(3600*24) #convert from days to seconds

print 'calculating temporal averages'
#2D fields
CRH_tmp = CRH.reshape(nt2D/ntave2D, ntave2D, nx, ny)
CRH_tave = np.mean(CRH_tmp, axis=1)

#3D fields (6-hourly output)
#reshape to take daily averages
qv_tmp = qv.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
T_tmp = T.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
w_tmp = w.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
rho_tmp = rhoz.reshape(nt3D/ntave3D, ntave3D, nz)
Qr_tmp = Qr.reshape(nt3D/ntave3D, ntave3D, nz)

qv_tave = np.mean(qv_tmp, axis=1)
T_tave = np.mean(T_tmp, axis=1)
w_tave = np.mean(w_tmp, axis=1)
Qr_tave = np.mean(w_tmp, axis=1)
rhoz_tave = np.mean(rho_tmp, axis=1)

z3D = np.ones((x.size, y.size, z.size))
z3D[:,:,:]=z
z3D = np.transpose(z3D)

#width of block in units of dx
#dx = 3 km, so db = 16*3 = 48 km
db=16

#look at last 40 days of simulation (only have 3D data from day 90 to day 130)
times=np.arange(-130,0)
nblocks = (nx/db)*(ny/db)

CRHs = np.zeros((times.size, nblocks))
whatsort = np.zeros(CRHs.shape)
Qrhatsort = np.zeros(CRHs.shape)

for i, ti in enumerate(times):
    
    #MSE and FMSE calculations
    print 'calculating vertically integrated velocity'
    print ti
    what = vertint(w_tave[ti,:,:,:], p)
    Qrhat = vertint(Qr_tave[ti,:,:,:], p)
    CRH_whatsort = blocksort2D(CRH_tave[ti,:,:], what, db)
    CRH_Qrhatsort = blocksort2D(CRH_tave[ti,:,:], Qrhat, db)
    CRHs[i,:] = CRH_whatsort.keys()
    whatsort[i,:] = CRH_whatsort.values()
    Qrhatsort[i,:] = CRH_Qrhatsort.values()
    
CRHranks = np.arange(np.size(CRHs[0,:]))
rankss, tt = np.meshgrid(CRHranks, times)
#convert tt to days
tt = tt+t3D[-1]

t = -1
a=1.5
CRHcrit = np.mean(CRHs[t, :]) + a*np.std(CRHs[t, :])
moist_indices = CRHs[t,:] > CRHcrit
dry_indices = CRHs[t,:] < CRHcrit

whatbar_moist = np.mean(whatsort[t, moist_indices])
whatbar_dry = np.mean(whatsort[t, dry_indices])

Qrhatbar_moist = np.mean(Qrhatsort[t, moist_indices])
Qrhatbar_dry = np.mean(Qrhatsort[t, dry_indices])

print 'domain size {:3.0f} km^2, day {:3.0f}'.format(domsize, t3D[-1])
print 'horizontal mean of density-weighted vertical integral of vertical velocity, moist region = {:3.2f} kg m^-2 s^-1'.format(whatbar_moist)
print 'horizontal mean of density-weighted vertical integral of vertical velocity, dry region = {:3.2f} kg m^-2 s^-1'.format(whatbar_dry)
print 'horizontal mean of density-weighted vertical integral of net radiative heating, moist region = {:3.2f} J m^-2 day^-1'.format(Qrhatbar_moist)
print 'horizontal mean of density-weighted vertical integral of net radiative heating, dry region = {:3.2f} J m^-2 day^-1'.format(Qrhatbar_dry)

plt.figure(1)
v = np.arange(-400, 400, 25)
cs = plt.contourf(rankss, tt, whatsort, v, cmap=cm.RdBu_r, extend='both')
cs.cmap.set_over('k')
plt.xlabel('CRH rank')
plt.ylabel('time (days)')
plt.title('$\hat w$ (kg m$^{-2}$ s$^{-1}$)')
plt.colorbar(cs)
plt.savefig(fout + 'CRHsort_whatsort_day{0}to{1}.pdf'.format(tt[0, 0], tt[-1,0]))
plt.close()

plt.figure(2)
cs = plt.contourf(rankss, tt, whatsort, v, cmap=cm.RdBu_r, extend='both')
cs.cmap.set_over('k')
plt.xlabel('CRH rank')
plt.ylabel('time (days)')
plt.title('$\hat Q_r$ (J m$^{-2}$ day$^{-1}$)')
plt.colorbar(cs)
plt.savefig(fout + 'CRHsort_Qrhatsort_day{0}to{1}.pdf'.format(tt[0, 0], tt[-1,0]))
plt.close()




