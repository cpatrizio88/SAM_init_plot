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
fout = '/Users/cpatrizio/figures/SAM_RCE100days_96kmnudgeuv/'

nc_in = glob.glob(fpath + '*1500m*100days_nudgeuv.nc')[0]
nc_data = Dataset(nc_in)
nc_vars = nc_data.variables

z = nc_vars['z'][:]
p = nc_vars['p'][:]
t = nc_vars['time'][:]
dt = (t[-1] - t[-2]) #difference between OUT_STAT output in days

profile_names = ['RELH', 'PRECIP', 'THETAE', 'THETA', 'QCOND', 'QV', 'TABS', 'U', 'V', 'MSE', 'DSE', 'RADQR']
timeseries_names = ['PW', 'LHF', 'SHF', 'PREC', 'CAPE', 'LWNTOA', 'LWNT', 'SWNTOA', 'LWDS', 'SWNS', 'CLDLOW'] 

for i, profile_name in enumerate(profile_names):
    profile_var = nc_vars[profile_name]
    profile = profile_var[:]
    tt, zz = np.meshgrid(t, z)
    tt, pp = np.meshgrid(t, p)
    nave = int(24*1) #average the profiles over nave hours to get the daily time tendency
    tendaveprofile = np.sum(profile[-nave:,:], axis=0)/nave
    tendaveprofile_prev = np.sum(profile[-2*nave:-nave, :], axis=0)/nave
    ddtprofile = (tendaveprofile - tendaveprofile_prev)/(nave*dt)
    frac_ddtprofile = ddtprofile/tendaveprofile_prev
    
    #t-z plots
    plt.figure(1)
    ax = plt.gca()
    plt.contourf(tt, pp, np.transpose(profile), 30, cmap=cm.RdYlBu_r)
    ax.set_yscale('log')
    plt.yticks([1000, 500, 250, 100, 50, 20])
    ax.set_ylim(p[-1], p[0])
    ax.invert_yaxis()
    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.colorbar(label='{0} ({1})'.format(profile_name, profile_var.units), orientation='horizontal')
    plt.title('{0} ({1})'.format(profile_name, profile_var.units))
    plt.xlabel('t (days)')
    plt.ylabel('p (hPa)')
    plt.savefig(fout + 'profile{0}days_'.format(np.round(t[-1])) + profile_name + '_idealRCE.pdf')
    plt.clf()
    
    #difference between model start and model end profiles
    plt.figure(2)
    ax = plt.gca()
    index_mid = np.where(t >= 60)[0][0]
    t0profile = profile[0,:]
    tmidprofile = profile[index_mid, :]
    tendprofile = profile[-1,:]
    #get the profile from the previous day (assume time step in hours)
    tprevendprofile = profile[-nave,:]
    plt.plot(t0profile, p, 'k--', label='{0} at t = 0'.format(profile_name))
    plt.plot(tmidprofile, p, 'b-', label='{0} at t = {1} days'.format(profile_name, np.round(t[index_mid])))
    plt.plot(tendaveprofile_prev, p, 'r-', label='{:2.2f}-day time averaged {:s} at t = {:3.1f} days'.format(nave/24., profile_name, t[-nave]))
    plt.plot(tendaveprofile, p, 'k-', label='{:2.2f}-day time averaged {:s} at t = {:3.1f} days'.format(nave/24., profile_name, t[-1]))
    ax.set_yscale('log')
    plt.yticks([1000, 500, 250, 100, 50, 20])
    ax.set_ylim(p[-1], p[0])
    ax.invert_yaxis()
    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xlabel('{0} ({1})'.format(profile_name, profile_var.units))
    plt.ylabel('p (hPa)')
    plt.title('{0} ({1})'.format(profile_name, profile_var.units))
    plt.legend()
    plt.savefig(fout + 'evoprofile{0}days_'.format(np.round(t[-1]))  + profile_name + '_idealRCE.pdf')
    plt.clf()
    
    #plot the lapse rate
    if (profile_name == 'TABS'):
        gamma = np.zeros(z.shape)
        delz = np.diff(z)/1e3 #delta z in km
        delT = np.diff(profile, axis=1)
        gamma0 = delT[0,:]/delz
        gamma40days = delT[index_mid, :]
        gammaend = delT[-1,:]
        plt.figure()
        ax = plt.gca()
        plt.plot(gamma0, p[:-1], 'k--', label='lapse rate at t = 0'.format(profile_name))
        plt.plot(gamma40days, p[:-1], 'b-', label='lapse rate at t = {:3.1f} days'.format(t[index_mid]))
        plt.plot(gammaend, p[:-1], 'k-', label='lapse rate at t = {:3.1f} days'.format(t[-1]))
        ax.set_yscale('log')
        plt.yticks([1000, 500, 250, 100, 50, 20])
        ax.set_ylim(p[-1], p[0])
        ax.invert_yaxis()
        ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.xlabel('lapse rate (K/km)')
        plt.ylabel('p (hPa)')
        plt.legend()
        plt.savefig(fout + 'evoprofile{0}days_'.format(np.round(t[-1])) + 'gamma_idealRCE.pdf')
        plt.clf()

        
    #time tendency (in units of per day) at model end 
    #get the difference between time-averaged profiles near the model end
    plt.figure(3)
    ax = plt.gca()
    plt.plot(ddtprofile, p)
    ax.set_yscale('log')
    ax.set_ylim(p[-1], p[0])
    plt.yticks([1000, 500, 250, 100, 50, 20])
    ax.invert_yaxis()
    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xlabel('time tendency of {0} ({1} per day)'.format(profile_name, profile_var.units))
    plt.ylabel('p (hPa)')
    plt.title('time tendency of {:s} ({:s} per day) at end of model run (using {:2.2f}-day time averaged profiles at end of model run)'.format(profile_name, profile_var.units, nave/24.))
    plt.savefig(fout + 'ddt{0}days_'.format(np.round(t[-1])) + profile_name + '_idealRCE.pdf')
    plt.clf()
    
    #percent change tendency 
    plt.figure(4)
    ax = plt.gca()
    plt.plot(frac_ddtprofile*100., p)
    ax.set_yscale('log')
    plt.yticks([1000, 500, 250, 100, 50, 20])
    ax.set_ylim(p[-1], p[0])
    ax.invert_yaxis()
    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xlabel('percent change time tendency (% per day)')
    plt.ylabel('p (hPa)')
    plt.title('percent change time tendency of {:s} (% per day) at end of model run (using {:2.2f}-day time averaged profiles at end of model run)'.format(profile_name, nave/24.))
    plt.savefig(fout + 'percentddt{0}days_'.format(np.round(t[-1])) + profile_name + '_idealRCE.pdf')
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
    
    #plot net radiation at TOA
    if (ts_name == 'LWNTOA'):
        plt.figure(2)
        LWNTOA_ts = ts
        SWNTOA_ts = nc_vars['SWNTOA'][:]
        QNTOA_ts = SWNTOA_ts - LWNTOA_ts 
        plt.title('QNTOA (W/m2)')
        plt.xlabel('t (days)')
        plt.ylabel('QNTOA (W/m2)')
        plt.plot(t, QNTOA_ts, 'k')
        plt.savefig(fout + 'timeseries{0}days_'.format(np.round(t[-1])) + 'QNTOA' + '_idealRCE.pdf')
        plt.clf()
        


   
    





