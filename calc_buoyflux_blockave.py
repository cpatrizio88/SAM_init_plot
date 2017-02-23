from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
from thermolib.constants import constants
import gc
from SAM_init_plot.block_fns import blockave3D

c = constants()

foutSTAT = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggr150days_64vert_ubarzero_STAT/'
foutSTATvary = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STAT/'
#foutSTATvary = '/Users/cpatrizio/Google Drive/figures/'

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpathSTAT = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
fout = '/Users/cpatrizio/data/SST302/'


matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 16})

nc_inSTAT = glob.glob(fpathSTAT + '*512x512*3000m*140days*302K.nc')[0]
nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day90to130*302K.nc')[0]

nc_data3D = Dataset(nc_in3D)
nc_dataSTAT = Dataset(nc_inSTAT)
varis3D = nc_data3D.variables
varisSTAT = nc_dataSTAT.variables

t2D = varisSTAT['time'][:]
t3D = varisSTAT['time'][:]

nt3D = t3D.size
nt2D = t2D.size

aveperiod2D = 24
aveperiod3D = 4

p = varis3D['p'][:]
p = p*1e2

z = varis3D['z'][:]
x = varis3D['x'][:]
y = varis3D['y'][:]

nz = z.size
nx = x.size
ny = y.size

domsizes = [1536]

#start and end time in days (corresponding to first file, and last file for given domain size)
tinits = [91]
colors = ['k', 'r', 'g'] 

#tdouble = 90


buoyfluxbar = []

buoyflux = np.array([])

fname = glob.glob(fout + '{:d}km_buoyflux_*130QP.npy'.format(domsizes[0]))[0]

#for i, d in enumerate(domsizes):
#    print 'domain size', d
#    fnames = glob.glob(fout + '{:d}km_buoyflux*QP.npy'.format(d))
#    for f in fnames:
#        print 'loading', f
#        buoyflux_temp = np.load(f)
#        buoyflux_temp = buoyflux_temp.astype('float32')
#        buoyflux = np.vstack((buoyflux, buoyflux_temp)) if buoyflux.size else buoyflux_temp
#    gc.collect()

buoyflux = np.load(fname)
buoyflux = buoyflux.astype('float32')

ntimes = buoyflux.shape[0]
delz = np.diff(z)
delz2D = np.zeros((ntimes, nz-1))
delz2D[:,:] = delz
buoyfluxbar_zt = np.mean(np.mean(buoyflux, axis=2), axis=2)
buoyfluxbar_int = np.sum(np.multiply(delz, buoyfluxbar_zt), axis=1)
buoyfluxbar = np.append(buoyfluxbar, buoyfluxbar_int)
times = np.arange(tinits[0], tinits[0]+ntimes)
        
#W = varis3D['W'][:]
#W = W.reshape(nt3D/aveperiod3D, aveperiod3D, nz, nx, ny)

#W_tave = np.mean(W, axis=1)

#W_anom = np.zeros(W_tave.shape)

W_prime = np.zeros((buoyfluxbar.shape[0], nz-1, nx, ny))
thetav_prime = np.zeros((buoyfluxbar.shape[0], nz-1, nx, ny))

db=16

buoyfluxblockbar = np.zeros(buoyfluxbar.shape)

print 'calculating block-averaged buoyancy flux'

for ti, t in enumerate(times):
    i = ti+1
    print 'day', t
    
    t2 = t*aveperiod2D
    t3 = i*aveperiod3D
    W = varis3D['W'][t3-aveperiod3D:t3,:,:,:]
    W_tave = np.mean(W, axis=0)
    
    theta = varisSTAT['THETA'][t2-aveperiod2D:t2,:]
    theta_bar = np.mean(theta, axis=0)
    
    T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
    rhoz = p/(c.Rd*np.mean(np.mean(T, axis=3), axis=2))
    rhoz_tave = np.mean(rhoz, axis=0)
    
    fac = (rhoz_tave*c.g/theta_bar)
    
    for k, zlev in enumerate(z[:-1]):
        #W_zbar = np.mean(np.mean(W_tave[k,:,:], axis=0), axis=0)
        #W_prime[ti,k,:,:] = W_tave[k,:,:] - W_zbar
        thetav_prime[ti,k,:,:] = buoyflux[ti,k,:,:]/(fac[k]*W_tave[k,:,:])
    
    #W_primeblock = blockave3D(W_prime[ti,:,:,:], db)
    W_block = blockave3D(W_tave, db)
    thetav_primeblock = blockave3D(thetav_prime[ti,:,:,:], db)
    
    W_block = W_block[:-1,:,:]
    
    #fac = (rhoz_tave*c.g/theta_bar)
    
    #fac3D = np.zeros((W_primeblock.shape)).T
    fac3D = np.zeros((W_block.shape)).T
    fac3D[:,:,:] = fac[:-1]
    fac3D = fac3D.T
        
    #buoyfluxblock = np.multiply(W_primeblock, thetav_primeblock)
    buoyfluxblock = np.multiply(W_block, thetav_primeblock)
    buoyfluxblock = np.multiply(buoyfluxblock, fac3D)
    buoyfluxblockbar_zt = np.mean(np.mean(buoyfluxblock, axis=2), axis=1)

    temp = np.multiply(delz, buoyfluxblockbar_zt)
    #TODO: WHERE ARE INFINITE/NAN VALUES COMING FROM? SURFACE NAN COMING FROM W = 0.
    buoyfluxblockbar[ti] = np.sum(temp[np.isfinite(temp)], axis=0)
    
    print 'buoyfluxblock mean', buoyfluxblockbar[ti]
    print 'buoyflux mean', buoyfluxbar[ti]
    
plt.figure(1)
plt.xlabel('time (days)')
plt.ylabel(r'vertically integrated domain-mean buoyancy flux (W/m$^2$)')
plt.title('domain-mean block-averaged buoyancy flux, domain width = {:d} km, block width = {:2.0f} km'.format(domsizes[0], (db*np.diff(x)[0])/1e3))
plt.plot(times, buoyfluxbar, 'kx-', label='mean')
plt.plot(times, buoyfluxblockbar, 'kx-', alpha=0.5, label='block-averaged mean')
plt.plot(times, buoyfluxbar - buoyfluxblockbar, 'bx-', alpha=0.5, label = 'sub-block scale mean')
plt.show()
    
    
    
        
#for each time, do 3D block averaging on W_prime and theta_vprime
# multiply together, and then do a domain average
# save block averaged wprime and thetavprime?
# save domain-averaged block-averaged buoyancy flux

#subtract domain-averaged block-averaged buoyancy flux from buoyfluxbar 
#to get the domain-averaged cumulus scale buoyancy flux



