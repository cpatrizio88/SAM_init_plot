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

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
fout = '/Users/cpatrizio/figures/SAM_RCE100days_nudgediff/'

nc_in1 = glob.glob(fpath + '*1500m*100days.nc')[0]
nc_in2 = glob.glob(fpath + '*5000m*100days.nc')[0]

delx1=1500
delx2=5000

nc_data1 = Dataset(nc_in1)
nc_data2 = Dataset(nc_in2)
nc_vars1 = nc_data1.variables
nc_vars2 = nc_data2.variables

z1= nc_vars1['z'][:]
z2 = nc_vars2['z'][:]



p1 = nc_vars1['p'][:]
p2 = nc_vars2['p'][:]

t1 = nc_vars1['time'][:]
t2 = nc_vars2['time'][:]

tt1, zz1 = np.meshgrid(t1, z1)
tt1, pp1 = np.meshgrid(t1, p1)

dt1 = (t1[-1] - t1[-2]) #difference between OUT_STAT output in days
dt2 = (t2[-1] - t2[-2]) 

tt2, zz2 = np.meshgrid(t2, z2)
tt2, pp2 = np.meshgrid(t2, p2)

profile_names = ['RELH', 'PRECIP', 'THETAE', 'THETA', 'QCOND', 'QV', 'TABS', 'U', 'V', 'MSE', 'DSE', 'RADQR']
timeseries_names = ['PW', 'LHF', 'SHF', 'PREC', 'CAPE', 'LWNTOA', 'LWNT', 'SWNTOA', 'LWDS', 'SWNS', 'CLDLOW'] 

nave = int(24*1) #average the profiles over nave hours to get the daily time tendency

for i, profile_name in enumerate(profile_names):
    profile_var1 = nc_vars1[profile_name]
    profile1 = profile_var1[:]

    profile_var2 = nc_vars2[profile_name]
    profile2 = profile_var2[:]
    
    tendaveprofile1 = np.sum(profile1[-nave:,:], axis=0)/nave
    tendaveprofile_prev1 = np.sum(profile1[-2*nave:-nave, :], axis=0)/nave
 
    tendaveprofile2 = np.sum(profile2[-nave:,:], axis=0)/nave
    tendaveprofile_prev2 = np.sum(profile2[-2*nave:-nave, :], axis=0)/nave
 
    #difference between model start and model end profiles
    plt.figure(1)
    ax = plt.gca()
    index_mid = np.where(t1 >= 60)[0][0]
    t0profile1 = profile1[0,:]
    tmidprofile1 = profile1[index_mid, :]
    tendprofile1 = profile1[-1,:]
    t0profile2 = profile2[0,:]
    tmidprofile2 = profile2[index_mid, :]
    tendprofile2 = profile2[-1,:]
    plt.plot(t0profile1, p1, 'k--', label='{0} at t = 0'.format(profile_name))
    #plt.plot(tmidprofile1, p1, 'b-', label='{0} at t = {1} days'.format(profile_name, np.round(t1[index_mid])))
    plt.plot(tendaveprofile_prev1, p1, 'r-', label='deltax = {:4.0f} m , {:2.2f}-day time averaged {:s} at t = {:3.1f} days'.format(delx1, nave/24., profile_name, t1[-nave]))
    plt.plot(tendaveprofile1, p1, 'k-', label='deltax = {:4.0f} m , {:2.2f}-day time averaged {:s} at t = {:3.1f} days'.format(delx1, nave/24., profile_name, t1[-1]))
    #plt.plot(t0profile2, p2, 'k--', alpha=0.2, label='ubar = 0, {0} at t = 0'.format(profile_name))
    #plt.plot(tmidprofile2, p2, 'b-', alpha=0.2, label='ubar = 0, {0} at t = {1} days'.format(profile_name, np.round(t1[index_mid])))
    plt.plot(tendaveprofile_prev2, p2, 'r-', alpha=0.2, label='deltax = {:4.0f} m, ubar = 0, {:2.2f}-day time averaged {:s} at t = {:3.1f} days'.format(delx2, nave/24., profile_name, t1[-nave]))
    plt.plot(tendaveprofile2, p2, 'k-', alpha=0.2, label='deltax = {:4.0f} m, ubar = 0, {:2.2f}-day time averaged {:s} at t = {:3.1f} days'.format(delx2, nave/24., profile_name, t1[-1]))
    ax.set_yscale('log')
    plt.yticks([1000, 500, 250, 100, 50, 20])
    ax.set_ylim(p2[-1], p2[0])
    ax.invert_yaxis()
    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xlabel('{0} ({1})'.format(profile_name, profile_var1.units))
    plt.ylabel('p (hPa)')
    plt.title('{0} ({1})'.format(profile_name, profile_var1.units))
    plt.legend()
    plt.savefig(fout + 'evoprofile_nudge_{0}days_'.format(np.round(t1[-1]))  + profile_name + '_idealRCE.pdf')
    plt.clf()
    

for i, ts_name in enumerate(timeseries_names):
    ts_var1 = nc_vars1[ts_name]
    ts1 = ts_var1[:]
    ts_var2 = nc_vars2[ts_name]
    ts2 = ts_var2[:]
    plt.figure(1)
    plt.plot(t1, ts1, '-k', label='deltax = {:4.0f} m'.format(delx1))
    plt.plot(t1, ts1, '-r', alpha=0.35, label='deltax = {:4.0f} m , ubar=0'.format(delx2))
    plt.title('{0} ({1})'.format(ts_name, ts_var1.units))
    plt.xlabel('t (days)')
    plt.ylabel('{0} ({1})'.format(ts_name, ts_var1.units))
    plt.title('{0} ({1})'.format(ts_name, ts_var1.units))
    plt.legend()
    plt.savefig(fout + 'timeseries_nudge_{0}days_'.format(np.round(t1[-1]))  + ts_name + '_idealRCE.pdf')
    plt.clf()
    
    #plot net radiation at TOA
    if (ts_name == 'LWNTOA'):
        plt.figure(2)
        LWNTOA_ts1 = ts1
        SWNTOA_ts1 = nc_vars1['SWNTOA'][:]
        QNTOA_ts1 = SWNTOA_ts1 - LWNTOA_ts1
        LWNTOA_ts2 = ts2
        SWNTOA_ts2 = nc_vars2['SWNTOA'][:]
        QNTOA_ts2 = SWNTOA_ts2 - LWNTOA_ts2
        plt.title('QNTOA (W/m2)')
        plt.xlabel('t (days)')
        plt.ylabel('QNTOA (W/m2)')
        plt.plot(t1, QNTOA_ts1, 'k', label='deltax = {:4.0f} m'.format(delx1))
        plt.plot(t2, QNTOA_ts2, 'r', alpha=0.35, label='deltax = {:4.0f} m , ubar=0'.format(delx2))
        plt.legend()
        plt.savefig(fout + 'timeseries_nudge_{0}days_'.format(np.round(t1[-1])) + 'QNTOA' + '_idealRCE.pdf')
        plt.clf()
    
    
    



