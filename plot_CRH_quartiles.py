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

domsize=768

nc_in3D = glob.glob(fpath3D + '*256x256*3000m*day90to130*302K.nc')[0]
nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

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

#averaging time period for 2D & 3D fields
ntave2D=24
ntave3D=4


#calculate CRH
PW = varis2D['PW'][:]
SWVP = varis2D['SWVP'][:]
CRH = PW/SWVP
CRH_tmp = CRH.reshape(nt2D/ntave2D, ntave2D, nx, ny)

#daily average
CRH_tave = np.mean(CRH_tmp, axis=1)

### plot mean q1, q2, q3, q4 profiles for given variable at a given time

#times in days
times = np.arange(int(t3D[0]), int(t3D[-1]))

########USER EDIT: variable and time of interest  
varname = 'QV'
vari = varis3D[varname]
print 'calculating daily average for variable {0}'.format(varname)
var_tmp = vari[:].reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
var_tave = np.mean(var_tmp, axis=1)

#time of interest
t = -15

#sort profiles by CRH
db=16
CRH_profsort = blocksort3D(CRH_tave[t,:,:], var_tave[t,:,:,:], db)
profsort = CRH_profsort[1]
CRHs_tavesort = CRH_profsort[0]

#calculating percentiles
CRHq1 = np.percentile(CRHs_tavesort, 25)
CRHq2 = np.percentile(CRHs_tavesort, 50)
CRHq3 = np.percentile(CRHs_tavesort, 75)
CRHq4 = np.percentile(CRHs_tavesort, 100)

profq1 = profsort[:, CRHs_tavesort < CRHq1]
profq2 = profsort[:, np.bitwise_and(CRHs_tavesort > CRHq1, CRHs_tavesort < CRHq2)]
profq3 = profsort[:, np.bitwise_and(CRHs_tavesort > CRHq2, CRHs_tavesort < CRHq3)]
profq4 = profsort[:, CRHs_tavesort > CRHq3]

profq1bar = np.mean(profq1, axis=1)
profq2bar = np.mean(profq2, axis=1)
profq3bar = np.mean(profq3, axis=1)
profq4bar = np.mean(profq4, axis=1)

plt.figure(1)
ax=plt.gcf().gca()
ax.set_yscale('log')
plt.yticks([1000, 500, 250, 100, 50, 20])
ax.set_ylim(p[-1]/1e2, p[0]/1e2)
ax.invert_yaxis()
ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.plot(profq1bar, p/1e2, label='q1')
plt.plot(profq2bar, p/1e2, label='q2')
plt.plot(profq3bar, p/1e2, label='q3')
plt.plot(profq4bar, p/1e2, label='q4')
plt.legend()
plt.xlabel('{:s} ({:s})'.format(varname, vari.units.strip()))
plt.ylabel('p (hPa)')
plt.title(r'daily-averaged {:s} profiles sorted by CRH quartiles, domain size {:3.0f} km$^2$'.format(varname, domsize))
plt.savefig(fout+ 'CRHsort_{:s}quartiles.pdf'.format(varname))
plt.close()