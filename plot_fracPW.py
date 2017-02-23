from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker
import matplotlib
from SAM_init_plot.block_fns import blockave3D, blockave2D
from thermolib.constants import constants

c = constants()

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 2})

plt.style.use('seaborn-white')
fpathSTAT = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
fpath2D =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fpath3D = '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
foutdata = '/Users/cpatrizio/data/SST302/'
#fout =  '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_MAPSNEW/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to140_64vert_ubarzero_MAPSNEW/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to150_64vert_ubarzero_MAPSNEW/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STAT/'


#nc_STAT = glob.glob(fpathSTAT + '*256x256*3000m*250days*302K.nc')[0]
#nc_in = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
 
nc_STAT = glob.glob(fpathSTAT + '*512x512*3000m*195days*302K.nc')[0]
nc_in = glob.glob(fpath2D + '*512x512*3000m*day090to140*302K.nc')[0]
nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day090to130*302K.nc')[0]

#nc_STAT = glob.glob(fpathSTAT + '*1024x1024*3000m*170days*302K.nc')[0]
#nc_in = glob.glob(fpath2D + '*1024x1024*3000m*day090to130*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]

#domsize=768
domsize=1536
#domsize=3072

nave2D = 24
nave3D = 4

nc_data= Dataset(nc_in)
nc_dataSTAT = Dataset(nc_STAT)
nc_data3D = Dataset(nc_in3D)
varis2D = nc_data.variables
varis3D = nc_data3D.variables
varisSTAT = nc_dataSTAT.variables
times2D = varis2D['time'][:]
times3D = varis3D['time'][:]

x=varis2D['x'][:]
y=varis2D['y'][:]

z = varis3D['z'][:]
p = varis3D['p'][:]
p=p*1e2

nx=x.size
ny=y.size
nz=z.size


#rho = p*c.Rd*varisSTAT['TABS'][:]

theta = varisSTAT['THETA'][:]

QV = varis3D['QV'][:]
QV = QV*1e-3

PW_tropbar = []
PW_BLbar = []
times = []
#z_BL = []

BLi = np.where(z > 1500)[0][0]
#BLi = 0

nc_3Dfs = glob.glob(fpath3D + '*512x512*.nc')


for nc_f in nc_3Dfs:
    print 'loading', nc_f
    nc_data3D = Dataset(nc_f)
    varis3D = nc_data3D.variables
    times3D = varis3D['time'][:]
    #times_temp = np.arange(times3D[0], times3D[-1], )
    times_temp = times3D
    
    QV = varis3D['QV'][:]
    QV = QV*1e-3

    
    PW_tropbar_temp = np.zeros(times_temp.size)
    PW_BLbar_temp = np.zeros(times_temp.size)
    z_BL_temp = np.zeros(times_temp.size)
    
    for i, t in enumerate(times_temp):
        #t2 = (i+1)*nave2D
        #t3 = (i+1)*nave3D
        #theta_tave = np.mean(theta[t2-nave2D:t2,:], axis=0)
        #theta_tave = theta[i,:]
        #dthetadz_crit = 2e-3
        #dthetadz = (theta_tave[1:]-theta_tave[:-1])/np.diff(z)
        #BLi = np.where(dthetadz > dthetadz_crit)[0][0]
        
        #z_BL_temp[i] = z[BLi]
        
        #QV_tave = np.mean(QV[t3-nave3D:t3,:,:,:], axis=0)
        QV_tave = QV[i,:,:,:]
        
        diffp3D = np.zeros((ny, nx, nz-1))
        diffp3D[:,:,:] = -1*np.diff(p)
        diffp3D=diffp3D.T
        
        PW_BL = 1./(c.rhol*c.g)*np.sum(np.multiply(QV_tave[:BLi,:,:], diffp3D[:BLi,:,:]), axis=0)
        PW_trop = 1./(c.rhol*c.g)*np.sum(np.multiply(QV_tave[BLi:-1,:,:], diffp3D[BLi:,:,:]), axis=0)
        
        PW_tropbar_temp[i] = np.mean(PW_trop)
        PW_BLbar_temp[i] = np.mean(PW_BL)
    
    times = np.append(times, times_temp)
    PW_tropbar = np.append(PW_tropbar, PW_tropbar_temp)
    PW_BLbar = np.append(PW_BLbar, PW_BLbar_temp)
    #z_BL = np.append(z_BL, z_BL_temp)
        
PW_tropbar = PW_tropbar*1e3
PW_BLbar = PW_BLbar*1e3

PWstat = varisSTAT['PW'][:]
timesSTAT = varisSTAT['time'][:]

plt.figure(1)
plt.plot(times, PW_BLbar, 'k-.', alpha=0.7, label='PW below z = {:3.0f} m'.format(z[BLi]))
plt.plot(times, PW_tropbar, 'k--', alpha=0.7, label='PW free troposphere')
plt.plot(times, PW_BLbar + PW_tropbar, color='k', label = 'total, calculate')
plt.plot(timesSTAT, PWstat, color = 'r', label= 'total, from STAT')
plt.xlabel('time (days)')
plt.ylabel('PW (mm)')
plt.title('daily-mean, domain-averaged PW')
plt.legend()
plt.savefig(fout + '{:d}km_fracPW_{:3.0f}days.pdf'.format(domsize, times[-1]))
plt.close()


#plt.figure(2)
#plt.plot(times, z_BL)
#plt.xlabel('time (days)')
#plt.ylabel('BL height (m)')
#plt.title(r'boundary layer height')
#plt.legend()
#plt.show()
    
    
    

