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

#fout = '/Users/cparizio/Google Drive/figures/SST302/vary

#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to130_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to150_64vert_ubarzero_MOISTDRYPROFS/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_MOISTPROFS/'


#nc_in = glob.glob(fpath + '*256x256*3000m*day230to250*302K.nc')[0]
#nc_in2 = glob.glob(fpath + '*512x512*3000m*day180to195*302K.nc')[0]
#nc_in3 = glob.glob(fpath + '*1024x1024*3000m*day170to180*302K.nc')[0]   

domsize=768

ntave2D=24
ntave3D=4


if domsize == 768:
    nave=10
elif domsize == 1536:
    nave=5
else:
    nave=3


#nc_ins = glob.glob(fpath + '*256x256*3000m*day230to250*302K.nc')
#nc_ins = glob.glob(fpath + '*512x512*3000m*302K.nc')
#nc_ins = glob.glob(fpath + '*512x512*3000m*180*302K.nc')
#nc_ins = glob.glob(fpath + '*1024x1024*3000m*180*302K.nc')
nc_ins = glob.glob(fpath + '*256x256*3000m*230*302K.nc')

sigma = np.array([np.empty(64)])
sigma = np.array([])
times = np.array([])

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
    
    sigmas = np.zeros((t3D.size, z.size))
    
    for i, t in enumerate(t3D):
        print t
        dely=np.diff(y)[0]
        W = varis3D['W'][i,:,:,:]
        #T = varis3D['TABS'][i,:,:,:]
        #Wtave = np.mean(W, axis=0)
        
        #EDIT: this critical w threshold determines where convective region is. 
        #What is a good way to choose this value? (want to neglect gravity wave perturbations)
        totpoints=nx*ny
        
        #rhoz = p/(c.Rd*np.mean(np.mean(T, axis=2), axis=1))
        #rhoz_tave = np.mean(rhoz, axis=0)
        
        sigma_temp = np.zeros(z.size)
    
        db=16
        Wcrit=0.01
        
        #block averaging
        Wblock = blockave3D(W, db)
        nxblock = nx // db
        nyblock = ny // db
        totpoints = nxblock*nyblock
        
        
        for j, zlev in enumerate(z):
            Wz = Wblock[j,:,:]
            sigma_temp[j] = len(Wz[Wz >= Wcrit])/(1.*totpoints)
        
        sigmas[i,:] = sigma_temp
    
    times = np.concatenate((times, t3D))
    sigma = np.vstack((sigma, sigmas)) if sigma.size else sigmas

        
zz, tt = np.meshgrid(z, times)
#sigma = sigma[1:,:]

#vals = np.arange(0, 0.09, 1e-4)
#vals = np.arange(0, 0.3, 0.005)
vmin=0
vmax=0.20
#vmax=0.3#
#vmax =0.05


t=-1
aveperiod3D = nave*ntave3D

sigma_end = np.mean(sigma[t-aveperiod3D:t,:],axis=0)
z = varis3D['z'][:]
    
plt.figure(1)
ax = plt.gcf().gca()
cs = plt.pcolormesh(tt, zz/1e3, sigma, vmin=vmin, vmax=vmax, cmap=cm.RdYlBu_r)
ax.set_xlim([tt[0,0], tt[-1,0]])
ax.set_ylim((0, 20))
plt.xlabel('time (days)')
plt.ylabel('height (km)')
plt.title('domain width = {:d} km'.format(domsize))
cb = plt.colorbar()
if db == 1:
    plt.title(r'convective fractional area, $w_{{c}}$ = {:2.3f} m/s'.format(Wcrit))
    ax.set_xlabel(r'$\sigma$')
else:
    plt.title(r'convective fractional area, $w_{{c}}$ = {:2.3f} m/s, block-averaging over ({:2.0f} km)$^2$'.format(Wcrit, db*(np.diff(x)[0])/1e3))
if db == 1:
    cb.set_label(r'$\sigma$')
else:
    cb.set_label(r'$\sigma_m$')
plt.savefig(fout + 'sigmavt_{:4.0f}km_{:3.0f}days_db{:d}_Wcrit{:3.3f}.pdf'.format(domsize, t3D[-1], db, Wcrit))
plt.close()

plt.figure(2)
ax = plt.gcf().gca()
if db == 1:
    plt.title(r'convective fractional area, $w_{{c}}$ = {:2.3f} m/s'.format(Wcrit))
    ax.set_xlabel(r'$\sigma$')
else:
    plt.title(r'convective fractional area, $w_{{c}}$ = {:2.3f} m/s, block-averaging over ({:2.0f} km)$^2$'.format(Wcrit, db*(np.diff(x)[0])/1e3))
    ax.set_xlabel(r'$\sigma_m$')
plt.plot(sigma_end, z/1e3, '-x', color='g', label='{:d} km, averaged over day {:2.0f} to {:2.0f}'.format(domsize, times[0], times[-1]))
#ax.set_xlabel(r'$\sigma$')
plt.legend(loc='best')
plt.ylabel('z (km)')
#axarr[0,].plot(sigma, z/1e3)
#axarr[0,].set_xlabel('sigma')
#axarr[0,].set_ylabel('z (km)')

ax.set_ylim(0, 20)
plt.savefig(fout + 'sigma_TEST_day250_{:d}day_db{:d}_Wc{:3.3f}.pdf'.format(nave, db, Wcrit))
plt.close()
     