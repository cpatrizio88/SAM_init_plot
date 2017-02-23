from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from thermolib.wsat import wsat
from thermolib.constants import constants
import SAM_init_plot.misc_fns 
from SAM_init_plot.block_fns import blockave3D
from SAM_init_plot.misc_fns import fracclusterarea

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (20, 10)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 16})

c = constants()

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'

#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to130_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to150_64vert_ubarzero_MOISTDRYPROFS/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_MOISTPROFS/'


nc_in = glob.glob(fpath + '*256x256*3000m*day230to250*302K.nc')[0]
#nc_in2 = glob.glob(fpath + '*512x512*3000m*day180to195*302K.nc')[0]
#nc_in3 = glob.glob(fpath + '*1024x1024*3000m*day170to180*302K.nc')[0]   

domsize=768

ntave2D=24
ntave3D=4


if domsize == 768:
    nave=20
elif domsize == 1536:
    nave=15
else:
    nave=5


nc_ins = glob.glob(fpath + '*256x256*3000m*302K.nc')[0]

sigma = []
times = []

for nc_f in nc_ins:

    nc_data = Dataset(nc_f)
    varis3D = nc_data.variables
    t3D = varis3D['time'][:]
    z = varis3D['z'][:]
    x = varis3D['x'][:]
    y = varis3D['y'][:]
    p = varis3D['p'][:]
    nx = x.size
    ny = y.size
    p = p*1e2
    
    for i, t in enumerate(t3D):
        dely=np.diff(y)[0]
        W = varis3D['W'][i,:,:,:]
        T = varis3D['TABS'][i,:,:,:]
        #Wtave = np.mean(W, axis=0)
        
        #EDIT: this critical w threshold determines where convective region is. 
        #What is a good way to choose this value? (want to neglect gravity wave perturbations)
        Wcrit=0.01
        totpoints=nx*ny
        
        rhoz = p/(c.Rd*np.mean(np.mean(T, axis=3), axis=2))
        rhoz_tave = np.mean(rhoz, axis=0)
        
        sigma_temp = np.zeros(z.size)
    
        db=16
        
        #block averaging
        Wblock = blockave3D(W, db)
        nxblock = nx // db
        nyblock = ny // db
        totpoints = nxblock*nyblock
        
        
        for i, zlev in enumerate(z):
            Wz = Wblock[i,:,:]
            sigma_temp[i] = len(Wz[Wz >= Wcrit])/(1.*totpoints)
        
        sigma = np.hstack(sigma, sigma_temp)
    
    times = np.concatenate(times, t3D)
    
tt, zz = np.meshgrid(times, z.size)
    
plt.figure(1)
plt.pcolormesh(tt, zz/1e3, sigma, cmap=cm.viridis)
plt.xlabel('time (days)')
plt.ylabel('height (km)')
cb = plt.colorbar()
cb.setlabel(r'$\sigma$')
plt.show()
            
     