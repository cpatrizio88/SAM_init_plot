from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker
import matplotlib

matplotlib.rcParams.update({'font.size': 24})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')


fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_STAT/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggr130days_64vert_ubarzero_STAT/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggr150days_64vert_ubarzero_STAT/'

domsize=768
#domsize=1536
#domsize = 3072

nc_in = glob.glob(fpath + '*256x256*3000m*250days_302K.nc')[0]
#nc_in = glob.glob(fpath + '*512x512*3000m*195days*302K*.nc')[0]
#nc_in = glob.glob(fpath + '*1024x1024*3000m*190days*302K*.nc')[0]


nc_data = Dataset(nc_in)
nc_vars = nc_data.variables

z = nc_vars['z'][:]
p = nc_vars['p'][:]
t = nc_vars['time'][:]
dt = (t[-1] - t[-2]) #difference between OUT_STAT output in days
tdouble = 90
#tsmallend=250

nstat=24

profile_names = ['RELH', 'PRECIP', 'THETAE', 'THETA', 'QN', 'QP', 'QCL', 'QI', 'QPL', 'QPI', 'QV', 'TABS', 'U', 'V', 'MSE', 'DSE', 'RADQR']
timeseries_names = ['PW', 'LHF', 'SHF', 'PREC', 'CAPE', 'LWNTOA', 'LWNT', 'SWNTOA', 'LWNS', 'SWNS', 'CLDLOW', 'CLDHI'] 

profile_names = ['QV', 'RADQR', 'QN', 'RELH', 'DSDZ']
#profile_names = ['DSDZ']

for i, profile_name in enumerate(profile_names):
    if profile_name == 'DSDZ':
        s = nc_vars['DSE']
        diffz = np.diff(z)
        diffs = np.diff(s, axis=1)
        diffz2D = np.zeros(diffs.shape)
        diffz2D[:,:]=diffz
        DSDZ = diffs/diffz
        profile = DSDZ
        units = 'J/m'
        z = z[:profile.shape[1]]
        p = p[:profile.shape[1]]
    else:
        profile_var = nc_vars[profile_name]
        profile = profile_var[:]
        units = profile_var.units
    
    if profile_name == 'DSDZ':
        titlename = r'$\frac{ds}{dz}$'
    elif profile_name == 'QV':
        titlename = r'$q_v$'
    elif profile_name == 'RELH':
        titlename = 'RH'
    elif profile_name == 'RADQR':
        titlename = r'$Q_r$'
    elif profile_name == 'QN':
        titlename = r'$Q_n$'
    else:
        titlename = profile_name
        
    #profile_var = nc_vars[profile_name]
    #profile = profile_var[:]
    tt, zz = np.meshgrid(t, z)
    tt, pp = np.meshgrid(t, p)
    nave = int(24*10) #average the profiles over nave hours to get the daily time tendency
    tendaveprofile = np.sum(profile[-nave:,:], axis=0)/nave
    tendaveprofile_prev = np.sum(profile[-2*nave:-nave, :], axis=0)/nave
    ddtprofile = (tendaveprofile - tendaveprofile_prev)/(nave*dt)
    frac_ddtprofile = ddtprofile/tendaveprofile_prev
    
    tinitprofile = np.zeros(profile.shape)
    tdoubleprofile = np.zeros(profile.shape)
    tinitprofile[:,:] = np.mean(profile[:nave,:], axis=0)
    tdoubleprofile[:,:] = np.mean(profile[tdouble*nstat-nave:tdouble*nstat+nstat,:], axis=0)
    
    tinitprofile[tinitprofile == 0] = np.nan 
    
    fracchange = (profile - tinitprofile)/tinitprofile
    fracchangedouble = (profile - tdoubleprofile)/tdoubleprofile
    
    fracchange[np.isnan(fracchange)] = 0
    fracchangedouble[np.isnan(fracchangedouble)] = 0

    
    
    
    #t-z plots
    plt.figure(1)
    ax = plt.gca()
    vmin=np.min(profile)
    vmax=np.max(profile)
    if vmin < 0:
       if np.abs(vmin) < np.abs(vmax):
          vmin = vmin - (vmax - np.abs(vmin))
       else:
          vmax = vmax + (np.abs(vmin) - vmax)
       cmap = cm.RdBu_r
    else:
       cmap = cm.GnBu
    plt.contourf(tt, zz/1e3, np.transpose(profile), 30, vmin=vmin, vmax=vmax, cmap=cmap)
    #ax.set_yscale('log')
    #plt.yticks([1000, 500, 250, 100, 50, 20])
    #ax.set_ylim(p[-1], p[0])
    #ax.invert_yaxis()
    #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.colorbar(label='{0} ({1})'.format(titlename, units), orientation='horizontal')
    plt.title(r'{:s} ({:s}), domain size ({:4.0f} km)$^2$'.format(titlename, units, domsize))
    plt.xlabel('t (days)')
    plt.ylabel('z (km)')
    plt.axvline(tdouble, color='b', alpha=0.8)
    plt.ylim((0, 16))
    plt.axvline(tdouble)
    plt.savefig(fout + 'profile{0}days_'.format(np.round(t[-1])) + profile_name + '_idealRCE.pdf')
    plt.clf()
    
    plt.figure(2)
    ax = plt.gca()
    vmin=np.min(fracchange)
    vmax=np.max(fracchange)
    if vmin < 0:
       if np.abs(vmin) < np.abs(vmax):
          vmin = vmin - (vmax - np.abs(vmin))
       else:
          vmax = vmax + (np.abs(vmin) - vmax)
       cmap = cm.RdBu_r
    else:
       cmap = cm.GnBu
    fracvals = np.linspace(-2,2,40)
    cs=plt.contourf(tt, zz/1e3, np.transpose(fracchange), fracvals, cmap=cmap, extend='both')
    #cs.cmap.set_over('k')
    #ax.set_yscale('log')
    #plt.yticks([1000, 500, 250, 100, 50, 20])
    #ax.set_ylim(p[-1], p[0])
    #ax.invert_yaxis()
    #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.colorbar(label='fractional change', orientation='horizontal')
    plt.title(r'fractional change in {:s}, domain size ({:4.0f} km)$^2$'.format(titlename, domsize))
    plt.xlabel('t (days)')
    plt.ylabel('z (km)')
    plt.axvline(tdouble, color='b', alpha=0.8)
    plt.ylim((0, 16))
    plt.axvline(tdouble)
    plt.savefig(fout + 'fracchangeprofile{0}days_'.format(np.round(t[-1])) + profile_name + '_idealRCE.png')
    plt.clf()
    
    doublefracvals = np.linspace(-1,1,30)
    
    plt.figure(7)
    cs=plt.contourf(tt, zz/1e3, np.transpose(fracchangedouble), doublefracvals, cmap=cmap, extend='both')
    plt.colorbar(label='fractional change', orientation='horizontal')
    plt.title(r'fractional change in {:s} after domain doubling, domain size ({:4.0f} km)$^2$'.format(titlename, domsize))
    plt.xlabel('t (days)')
    plt.ylabel('z (km)')
    plt.axvline(tdouble, color='b', alpha=0.8)
    plt.ylim((0, 16))
    plt.axvline(tdouble)
    plt.savefig(fout + 'tdoublefracchangeprofile{0}days_'.format(np.round(t[-1])) + profile_name + '_idealRCE.png')
    plt.clf()
    
    #difference between model start and model end profiles
    plt.figure(3)
    ax = plt.gca()
    #EDIT: TIME TO LOOK AT AVERAGE PROFILES
    ti = 169
    index_taggr = np.where(t >= ti)[0][0]
    t0 = nave*10
    t0 = 0
    if t0 == 0:
        t0profile = profile[0,:]
    else:
        t0profile = np.mean(profile[t0-nave:t0,:], axis=0)
    tmidprofile = np.mean(profile[index_taggr-nave:index_taggr,:], axis=0)
    tmidprofile_prev = np.mean(profile[index_taggr-2*nave:index_taggr-nave, :], axis=0)
    #tendprofile = profile[-1,:]
    #get the profile from the previous day (assume time step in hours)
    tprevendprofile = profile[-nave,:]
    #plt.plot(t0profile, z/1e3, 'k--', label='{0} at day {1}'.format(profile_name, np.round(t[t0])))
    plt.plot(tmidprofile_prev, z/1e3, 'k-', alpha=0.2, label='{0} at day {1}'.format(titlename, np.round(t[index_taggr-nave])))
    plt.plot(tmidprofile, z/1e3, 'k-', label='{0} at day {1}'.format(titlename, np.round(t[index_taggr])))
    #plt.plot(tendaveprofile_prev, p, 'r-', label='{:2.2f}-day time averaged {:s} at t = {:3.1f} days'.format(nave/24., profile_name, t[-nave]))
    #plt.plot(tendaveprofile, p, 'k-', label='{:2.2f}-day time averaged {:s} at t = {:3.1f} days'.format(nave/24., profile_name, t[-1]))
    #ax.set_yscale('log')
    #plt.yticks([1000, 500, 250, 100, 50, 20])
    #ax.set_ylim(p[-1], p[0])
    #ax.invert_yaxis()
    #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    if profile_name == 'QI':
       ax.set_xlim([0, 0.07])
    plt.xlabel('{0} ({1})'.format(titlename, units))
    plt.ylabel('z (km)')
    plt.ylim(0,16)
    plt.title(r'{:s} ({:s}), domain size ({:4.0f} km)$^2$'.format(titlename, units, domsize))
    plt.legend()
    plt.savefig(fout + 'evoprofile{0}days_'.format(np.round(t[index_taggr]))  + profile_name + '_idealRCE.pdf')
    plt.clf()
    
    
    
    #plot the lapse rate
    if (profile_name == 'TABS'):
        gamma = np.zeros(z.shape)
        delz = np.diff(z)/1e3 #delta z in km
        delT = np.diff(profile, axis=1)
        gamma0 = delT[0,:]/delz
        gamma40days = delT[index_taggr, :]/delz
        gammaend = delT[-1,:]/delz
        plt.figure()
        ax = plt.gca()
        plt.plot(gamma0, p[:-1], 'k--', label='lapse rate at t = 0'.format(profile_name))
        plt.plot(gamma40days, p[:-1], 'b-', label='lapse rate at t = {:3.1f} days'.format(t[index_taggr]))
        plt.plot(gammaend, p[:-1], 'k-', label='lapse rate at t = {:3.1f} days'.format(t[-1]))
        #ax.set_yscale('log')
        #plt.yticks([1000, 500, 250, 100, 50, 20])
        #ax.set_ylim(p[-1], p[0])
        #ax.invert_yaxis()
        #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.xlabel('lapse rate (K/km)')
        plt.ylabel('p (hPa)')
        plt.ylim(0,16)
        plt.legend()
        plt.savefig(fout + 'evoprofile{0}days_'.format(np.round(t[index_taggr])) + 'gamma_idealRCE.pdf')
        plt.clf()

        
    ##time tendency (in units of per day) at model end 
    ##get the difference between time-averaged profiles near the model end
    #plt.figure(4)
    #ax = plt.gca()
    #plt.plot(ddtprofile, z/1e3)
    ##ax.set_yscale('log')
    ##ax.set_ylim(p[-1], p[0])
    ##plt.yticks([1000, 500, 250, 100, 50, 20])
    ##ax.invert_yaxis()
    ##ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #plt.xlabel('time tendency of {0} ({1} per day)'.format(titlename, units))
    #plt.ylabel('p (hPa)')
    #plt.ylim(0,16)
    #plt.title(r'time tendency of {:s} ({:s} per day) at end of model run (using {:2.2f}-day time averaged profiles at end of model run), domain size ({:4.0f} km)$^2$ '.format(profile_name, profile_var.units, nave/24., domsize))
    #plt.savefig(fout + 'ddt{0}days_'.format(np.round(t[-1])) + profile_name + '_idealRCE.pdf')
    #plt.clf()
    
    #percent change tendency 
    #plt.figure(4)
    #ax = plt.gca()
    #plt.plot(frac_ddtprofile*100., p)
    #ax.set_yscale('log')
    #plt.yticks([1000, 500, 250, 100, 50, 20])
    #ax.set_ylim(p[-1], p[0])
    #ax.invert_yaxis()
    #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #plt.xlabel('percent change time tendency (% per day)')
    #plt.ylabel('p (hPa)')
    #plt.title('percent change time tendency of {:s} (% per day) at end of model run (using {:2.2f}-day time averaged profiles at end of model run)'.format(profile_name, nave/24.))
    #plt.savefig(fout + 'percentddt{0}days_'.format(np.round(t[-1])) + profile_name + '_idealRCE.pdf')
    #plt.clf()
    
for i, ts_name in enumerate(timeseries_names):
    
    
    if ts_name == 'PREC':
        titlename='P'
    elif ts_name == 'LWNTOA':
        titlename=r'LW$_{net,TOA}$'
    elif ts_name == 'SWNTOA':
        titlename='SW$_{net,TOA}$'
    elif ts_name == 'CLDLOW':
        titlename = 'Low cloud fraction'
    elif ts_name == 'CLDHI':
        titlename = 'High cloud fraction'
    else:
        titlename = ts_name
        
    ts_var = nc_vars[ts_name]
    ts = ts_var[:]
    plt.figure(1)
    plt.plot(t, ts, '-k')
    plt.title(r'{:s} ({:s}), domain size ({:4.0f} km)$^2$'.format(titlename, ts_var.units, domsize))
    plt.xlabel('time (days)')
    plt.ylabel('{0} ({1})'.format(titlename, ts_var.units))
    plt.axvline(tdouble)
    plt.savefig(fout + 'timeseries{0}days_'.format(np.round(t[-1]))  + ts_name + '_idealRCE.pdf')
    plt.clf()
    
    #plot net radiation at TOA
    if (ts_name == 'LWNTOA'):
        plt.figure(2)
        LWNTOA_ts = ts
        SWNTOA_ts = nc_vars['SWNTOA'][:]
        QNTOA_ts = SWNTOA_ts - LWNTOA_ts 
        plt.title(r'QNTOA (W/m2), domain size ({:4.0f} km)$^2$'.format(domsize))
        plt.xlabel('time (days)')
        plt.ylabel('QNTOA (W/m2)')
        plt.plot(t, QNTOA_ts, 'k')
        plt.axvline(tdouble)
        plt.savefig(fout + 'timeseries{0}days_'.format(np.round(t[-1])) + 'QNTOA' + '_idealRCE.pdf')
        plt.clf()
        


   
    





