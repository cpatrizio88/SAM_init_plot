from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker
import matplotlib

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath2D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath1D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/figures/SST302/SAM_aggr100days_12288km_64vert_ubarzero_XZ2D/'

nc_in1D = glob.glob(fpath1D + '*4096x64*3000m*100days*302K.nc')[0]

nc_in2D = glob.glob(fpath2D + '*4096x64*3000m*100days*302K.nc')[0]

nc_data2D = Dataset(nc_in2D)
nc_data1D = Dataset(nc_in1D)
varis1D = nc_data1D.variables
varis2D = nc_data2D.variables
t2D = varis2D['time'][:]

x = varis2D['x'][:]
z = varis2D['z'][:]
p = varis2D['p'][:]

PW = varis1D['PW'][:]
u = varis2D['U'][:]
w = varis2D['W'][:]

xx, zz = np.meshgrid(x, z)
#xx, pp = np.meshgrid(x, p)

#average fields over some time period
nave=8*1 #time step of 2D fields is 3 hours (8 steps in a day for 2D simulations)
delt = np.diff(t2D)[0]*24 #time slice in hours

t=8*55

u_tave = np.mean(u[t-nave:t,:,:],axis=0)
w = varis2D['W'][:]
w_tave = np.mean(w[t-nave:t,:,:], axis=0)

field_names = ['TABS', 'QRAD' ,'U', 'W', 'QN', 'QP', 'QV']

for name in field_names:
    plt.figure()
    ax = plt.gca()
    vari = varis2D[name]
    avefieldt = np.mean(vari[t-nave:t,:,:],axis=0)
    fieldbar = np.mean(avefieldt)
    anomavefieldt = avefieldt - fieldbar
    if name == 'U':
        plt.contour(xx/1e3, zz/1e3, avefieldt, 6, colors='k', alpha=0.3)
        plt.contourf(xx/1e3, zz/1e3, avefieldt, 15, cmap = cm.RdYlBu_r)
    #TODO: work on quiver plotting..
    #elif name == 'QRAD':
    #    q = plt.quiver(xx[::4, ::64]/1e3, zz[::4, ::64]/1e3, u_tave[::4,::64], w_tave[::4,::64], scale=100, alpha=0.8, headwidth=2, angles='xy', zorder=1)
    #    plt.contourf(xx/1e3, zz/1e3, avefield_t, 30, cmap = cm.RdYlBu_r, zorder=0)
    #    plt.quiverkey(q, np.min(x)/1e3+30, np.max(t), 10, "10 m/s",coordinates='data',color='k')
    #elif name == 'QN':
    #    plt.contourf(xx/1e3, zz/1e3, avefield_t, np.linspace(0, 0.01, 65), cmap = cm.RdYlBu_r)
    else:
        plt.contourf(xx/1e3, zz/1e3, avefieldt, 15, cmap = cm.RdYlBu_r)
    plt.xlabel('x (km)')
    plt.ylabel('z (km)')
    #plt.ylabel('p (hPa)')
    #ax.set_yscale('log')
    #plt.yticks([1000, 950, 500, 250, 150, 50, 20])
    #ax.set_ylim(p_trop*1e-2, p[0])
    #ax.invert_yaxis()
    #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.title(r'{:s} [{:s}] temporally-averaged over days {:2.0f} to {:2.0f}, $\overbar${:s} = {:3.1f} {:s}'.format(name, vari.units.strip(), t2D[t-nave], t2D[t], name, fieldbar, vari.units.strip()))
    plt.grid()
    cb=plt.colorbar()
    cb.set_label('{:s}'.format(vari.units))
    plt.savefig(fout + '{:s}_day{:2.1f}to{:2.1f}.pdf'.format(name, t2D[t-nave], t2D[t]))
    plt.close()
    
#plot mean profiles for moist and dry region (vertical velocity profile, radiative warming rate profile, temperature profile, water vapor profile)

#plot boundary layer variables (relative humidity)

#plt.figure(2)
#plt.contour(xx[p > p_BL, :]/1e3, zz[p > p_BL, :]/1e3, qv_t_maxPW[p > p_BL, :], np.linspace(minqv, maxqv, 12), colors='k', alpha=0.5)
#plt.contourf(xx[p > p_BL, :]/1e3, zz[p > p_BL, :]/1e3, qv_t_maxPW[p > p_BL, :], np.linspace(minqv, maxqv, 12), cmap = cm.RdYlBu_r)
#plt.xlabel('x (km)')
#plt.ylabel('z (km)')
#plt.title('boundary layer water vapor mixing ratio (kg/kg) at y = y(max{var_x(PW)})')
#plt.colorbar()
#plt.show()

 








