from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib

matplotlib.rcParams.update({'font.size': 14})

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
fout = '/Users/cpatrizio/Dropbox/research/SAM figures/'

nc_in = glob.glob(fpath + '*.nc')[0]

nc_data = Dataset(nc_in)

nc_vars = nc_data.variables

z = nc_vars['z'][:]

p = nc_vars['p'][:]

t = nc_vars['time'][:]


profile_names = ['RELH', 'PRECIP', 'THETAE', 'THETA', 'QCOND', 'QV', 'TABS', 'U', 'V', 'MSE', 'DSE', 'RADQR']

timeseries_names = ['PW', 'LHF', 'SHF', 'PREC', 'CAPE', 'LWNTOA'] 


for i, profile_name in enumerate(profile_names):
    profile_var = nc_vars[profile_name]
    profile = profile_var[:]
    tt, zz = np.meshgrid(t, z)
    plt.figure(1)
    ax = plt.gca()
    plt.contourf(tt, zz/1000, np.transpose(profile), 20, cmap=cm.RdYlBu_r)
    plt.colorbar(label='{0} ({1})'.format(profile_name, profile_var.units), orientation='horizontal')
    plt.title('{0} ({1})'.format(profile_name, profile_var.units))
    plt.xlabel('t (days)')
    plt.ylabel('z (km)')
    plt.savefig(fout + profile_name + '_idealRCE{0}days.pdf'.format(np.round(t[-1])))
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
    plt.savefig(fout + ts_name + '_idealRCE{0}days.pdf'.format(np.round(t[-1])))
    plt.clf()
   
    





