from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib

matplotlib.rcParams.update({'font.size': 14})

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
fout = '/Users/cpatrizio/Dropbox/research/SAM RCE figures/'

nc_in = glob.glob(fpath + '*.nc')[0]

nc_data = Dataset(nc_in)

nc_vars = nc_data.variables

z = nc_vars['z'][:]

p = nc_vars['p'][:]

t = nc_vars['time'][:]

dt = t[-1] - t[-2]


profile_names = ['RELH', 'PRECIP', 'THETAE', 'THETA', 'QCOND', 'QV', 'TABS', 'U', 'V', 'MSE', 'DSE', 'RADQR']

timeseries_names = ['PW', 'LHF', 'SHF', 'PREC', 'CAPE', 'LWNTOA'] 


for i, profile_name in enumerate(profile_names):
    profile_var = nc_vars[profile_name]
    profile = profile_var[:]
    t0profile = profile[0,:]
    tendprofile = profile[-1, :]
    tt, zz = np.meshgrid(t, z)
    #t-z plots
    plt.figure(1)
    ax = plt.gca()
    plt.contourf(tt, zz/1000, np.transpose(profile), 20, cmap=cm.RdYlBu_r)
    plt.colorbar(label='{0} ({1})'.format(profile_name, profile_var.units), orientation='horizontal')
    plt.title('{0} ({1})'.format(profile_name, profile_var.units))
    plt.xlabel('t (days)')
    plt.ylabel('z (km)')
    plt.savefig(fout + 'profile{0}days_'.format(np.round(t[-1])) + profile_name + '_idealRCE.pdf')
    plt.clf()
    #difference between model start and model end profiles
    plt.figure(2)
    plt.plot(t0profile, z, label='{0} at t = 0'.format(profile_name))
    plt.plot(tendprofile, z, label='{0} at t = {1} days'.format(profile_name, np.round(t[-1])))
    plt.xlabel('{0} ({1})'.format(profile_name, profile_var.units))
    plt.ylabel('z (km)')
    plt.title('change in {0} ({1}) over {2} days'.format(profile_name, profile_var.units, np.round(t[-1])))
    plt.legend()
    plt.savefig(fout + 'diffprofile{0}days_'.format(np.round(t[-1]))  + profile_name + '_idealRCE.pdf')
    plt.clf()
    #time tendency (in units of per day) at model end 
    plt.figure(3)
    ddtprofile = (profile[-1, :] - profile[-2,:])/dt
    plt.plot(ddtprofile, z)
    plt.xlabel('time tendency of {0} ({1} per day)'.format(profile_name, profile_var.units))
    plt.ylabel('z (km)')
    plt.title('time tendency of {0} ({1} per day) at t = {2} days'.format(profile_name, profile_var.units, np.round(t[-1])))
    plt.savefig(fout + 'ddt{0}days_'.format(np.round(t[-1])) + profile_name + 'idealRCE.pdf')
    plt.clf()
  
    
for i, ts_name in enumerate(timeseries_names):
    ts_var = nc_vars[ts_name]
    ts = ts_var[:]
    plt.figure(1)
    plt.plot(t, ts, '-k')
    plt.title('{0} ({1})'.format(ts_name, ts_var.units))
    plt.xlabel('t (days)')
    plt.ylabel('{0} ({1})'.format(ts_name, ts_var.units))
    plt.title('{0} ({1})'.format(ts_name, ts_var.units))
    plt.savefig(fout + 'timeseries{0}days_'.format(np.round(t[-1]))  + ts_name + '_idealRCE.pdf')
    plt.clf()
    


   
    





