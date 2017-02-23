from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker
import matplotlib

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
times = varis['time'][:]

x = varis['x'][:]
y = varis['y'][:]
z = varis['z'][:]
p = varis['p'][:]

time=varis['time'][:]

xx, zz = np.meshgrid(x, z)

xxx, yyy = np.meshgrid(x, y)

#p = p*1e2 #convert from hPa to Pa
qv = varis['QV'][:]
qv = qv*1e-3 #convert from g/kg to kg/kg

p_BL = 950e2
p_trop = 150e2
p = np.concatenate(([p_s], p))

t=-1
qv_t = qv[t,:,:,:]

delp3D = np.ones((x.size, y.size, z.size))

delp = -np.diff(p)*1e2

delp3D[:,:,:] = delp

delp3D = np.transpose(delp3D)

#find the y coord where variance of PW is max.
#then plot different quantities for dry region/moist region (radiative heating, water vapor, vertical velocity, T)
PW = 1/(g*rho_w)*np.sum(np.multiply(delp3D, qv_t), axis=0)

minPW= np.min(PW)
maxPW= np.max(PW)

#plt.figure(1)
#plt.contour(xxx/1e3, yyy/1e3, PW, np.linspace(minPW, maxPW, 12), colors='k', alpha=0.5)
#plt.contourf(xxx/1e3, yyy/1e3, PW, np.linspace(minPW, maxPW, 12), cmap = cm.RdYlBu_r)
#plt.xlabel('x (km)')
#plt.ylabel('y (km)')
#plt.title('PW [m] at t = {:3.2f} days'.format(times[tplot]))
#plt.colorbar()
#plt.show()

#look at y-slices of different variables
#find the y-coord where the x-variance of PW is maximum
varPWx = np.var(PW, axis=1)
indx = np.where(varPWx == np.max(varPWx))[0]

#average fields over some time period
nave=4*10 #time step of 3D fields is 6 hours 
delt = np.diff(times)[0]*24 #time slice in hours

qv_t_yslice = np.squeeze(qv_t[:,:,indx])

aveqv_t = np.mean(qv[t-nave:t,:,:,:],axis=0)
aveqv_t_yslice = np.squeeze(aveqv_t[:,:,indx])

minqv = np.min(qv)
maxqv = np.max(qv)

plt.figure(1)
ax = plt.gca()
plt.contour(xx/1e3, zz/1e3, aveqv_t_yslice, np.linspace(minqv, maxqv, 24), colors='k', alpha=0.5)
plt.contourf(xx/1e3, zz/1e3, aveqv_t_yslice, np.linspace(minqv, maxqv, 24), cmap = cm.RdYlBu_r)
plt.xlabel('x (km)')
plt.ylabel('z (km)')
#plt.ylabel('p (hPa)')
#ax.set_yscale('log')
#plt.yticks([1000, 950, 500, 250, 150, 50, 20])
#ax.set_ylim(p_trop*1e-2, p[0])
#ax.invert_yaxis()
#ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.title('{:s} [{:s}] averaged over days {:2.1f} to {:2.1f} at y = y(max(var_x(PW)))'.format('QV', 'kg/kg', time[t-nave], time[t]))
plt.grid(zorder=2)
plt.colorbar()
plt.savefig(fout + 'QV_maxPWslice_day{:2.1f}to{:2.1f}.pdf'.format(times[t-nave], times[t]))
plt.close()

field_names = ['TABS', 'QRAD' ,'U', 'W', 'QN', 'QP']

for name in field_names:
   
    plt.figure()
    ax = plt.gca()
    vari = varis[name]
    avefield_t = np.mean(vari[t-nave:t,:,:,:],axis=0)
    avefield_t_yslice = np.squeeze(avefield_t[:,:,indx])
    #minfield = np.min(vari[:])
    #maxfield = np.max(vari[:])
    #plt.contour(xx/1e3, pp, field_t_yslice, 20, colors='k', alpha=0.5)
    if name == 'QN':
         plt.contourf(xx/1e3, zz/1e3, avefield_t_yslice, np.linspace(0, 0.01, 65), cmap = cm.RdYlBu_r)
    else:
         plt.contourf(xx/1e3, zz/1e3, avefield_t_yslice, 30, cmap = cm.RdYlBu_r)
    plt.xlabel('x (km)')
    plt.ylabel('z (km)')
    #plt.ylabel('p (hPa)')
    #ax.set_yscale('log')
    #plt.yticks([1000, 950, 500, 250, 150, 50, 20])
    #ax.set_ylim(p_trop*1e-2, p[0])
    #ax.invert_yaxis()
    #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.title(' {:s} [{:s}] temporally-averaged over days {:2.1f} to {:2.1f} at y = y(max(var_x(PW)))'.format(name, vari.units, time[t-nave], time[t]))
    plt.grid()
    plt.colorbar()
    plt.savefig(fout + '{:s}_maxPWslice_day{:2.1f}to{:2.1f}.pdf'.format(name, time[t-nave], time[t]))
    plt.close()
    
    #plot meridional average here
    
    plt.figure()
    ax = plt.gca()
    mavefield_t = np.mean(avefield_t, axis=2)
    if name == 'QN':
         plt.contourf(xx/1e3, zz/1e3, mavefield_t, np.linspace(0, 0.01, 65), cmap = cm.RdYlBu_r)
    else:
         plt.contourf(xx/1e3, zz/1e3, mavefield_t, 30, cmap = cm.RdYlBu_r)
    plt.xlabel('x (km)')
    plt.ylabel('z (km)')
    #plt.ylabel('p (hPa)')
    #ax.set_yscale('log')
    #plt.yticks([1000, 950, 500, 250, 150, 50, 20])
    #ax.set_ylim(p_trop*1e-2, p[0])
    #ax.invert_yaxis()
    #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.title('meridional average of {:s} [{:s}] temporally-averaged over days {:2.1f} to {:2.1f}'.format(name, vari.units, time[t-nave], time[t]))
    plt.grid()
    plt.colorbar()
    plt.savefig(fout + '{:s}_meridave_day{:2.1f}to{:2.1f}.pdf'.format(name, time[t-nave], time[t]))
    plt.close()
    
    

#plot mean profiles for moist and dry region (vertical velocity profile, radiative warming rate profile, temperature profile, water vapor profile)

#plot boundary layer variables (relative humidity)

#plot stream lines (Bretherton method of using CRH as horizontal coord vs. u,w ??) and overlay on RH





#plt.figure(2)
#plt.contour(xx[p > p_BL, :]/1e3, zz[p > p_BL, :]/1e3, qv_t_maxPW[p > p_BL, :], np.linspace(minqv, maxqv, 12), colors='k', alpha=0.5)
#plt.contourf(xx[p > p_BL, :]/1e3, zz[p > p_BL, :]/1e3, qv_t_maxPW[p > p_BL, :], np.linspace(minqv, maxqv, 12), cmap = cm.RdYlBu_r)
#plt.xlabel('x (km)')
#plt.ylabel('z (km)')
#plt.title('boundary layer water vapor mixing ratio (kg/kg) at y = y(max{var_x(PW)})')
#plt.colorbar()
#plt.show()

 








