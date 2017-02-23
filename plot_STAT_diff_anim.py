from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker
import matplotlib

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 2})

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STATANIM/'

nc_in1 = glob.glob(fpath + '*256x256*3000m*250days_302K.nc')[0]
nc_in2 = glob.glob(fpath + '*512x512*3000m*180days_302K.nc')[0]
nc_in3 = glob.glob(fpath + '*1024x1024*3000m*180days_302K.nc')[0]

delx1=3000
domsize1 = 768
delx2=3000
domsize2 = 1536
domsize3 = 3072

tdouble=90

nc_data1 = Dataset(nc_in1)
nc_data2 = Dataset(nc_in2)
nc_data3 = Dataset(nc_in3)
nc_vars1 = nc_data1.variables
nc_vars2 = nc_data2.variables
nc_vars3 = nc_data3.variables


z1= nc_vars1['z'][:]
z2 = nc_vars2['z'][:]
z3 = nc_vars3['z'][:]

p1 = nc_vars1['p'][:]
p2 = nc_vars2['p'][:]
p3 = nc_vars3['p'][:]

t1 = nc_vars1['time'][:]
t2 = nc_vars2['time'][:]
t3 = nc_vars3['time'][:]

tt1, zz1 = np.meshgrid(t1, z1)
tt1, pp1 = np.meshgrid(t1, p1)

dt1 = (t1[-1] - t1[-2]) #difference between OUT_STAT output in days
dt2 = (t2[-1] - t2[-2]) 
dt3 = (t3[-1] - t3[-2])

tt2, zz2 = np.meshgrid(t2, z2)
tt2, pp2 = np.meshgrid(t2, p2)

tt3, zz3 = np.meshgrid(t3, z3)
tt3, pp3 = np.meshgrid(t3, p3)

profile_names = ['RELH', 'PRECIP', 'THETAE', 'THETA', 'QN', 'QV', 'QCL', 'QP', 'QPI', 'QCOND', 'TABS', 'U', 'V', 'MSE', 'DSE', 'RADQR', ]
#timeseries_names = ['PW', 'LHF', 'SHF', 'PREC', 'CAPE', 'LWNTOA', 'LWNT', 'SWNTOA', 'LWDS', 'SWNS', 'CLDLOW', 'CLDHI'] 

profile_names = ['RELH']

nstat=24

nave = int(nstat) #average the profiles over nave hours to get the daily time tendency

times = np.arange(150, 180)

xmin=0
xmax=100

for ti, t in enumerate(times):
    t1i = t*nstat
    t2i = t*nstat
    t3i = t*nstat
    

    for i, profile_name in enumerate(profile_names):
        profile_var1 = nc_vars1[profile_name]
        profile1 = profile_var1[:]
    
        profile_var2 = nc_vars2[profile_name]
        profile2 = profile_var2[:]
        
        profile_var3 = nc_vars3[profile_name]
        profile3 = profile_var3[:]
        
        tendaveprofile1 = np.mean(profile1[t1i-nstat:t1i,:], axis=0)
        #tendaveprofile_prev1 = np.sum(profile1[-2*nave:-nave, :], axis=0)/nave
    
        tendaveprofile2 = np.mean(profile2[t2i-nstat:t2i,:], axis=0)
        #tendaveprofile_prev2 = np.sum(profile2[-2*nave:-nave, :], axis=0)/nave
        
        tendaveprofile3 = np.mean(profile3[t3i-nave:t3i,:], axis=0)
    
        #difference between model start and model end profiles
        plt.figure(1)
        ax = plt.gca()
        index_mid = np.where(t1 >= 60)[0][0]
        plt.plot(tendaveprofile1, z1/1e3, 'k-', label='{:3.0f} km'.format(domsize1))
        plt.plot(tendaveprofile2, z2/1e3, 'r-', label='{:4.0f} km'.format(domsize2))
        plt.plot(tendaveprofile3, z3/1e3, 'g-', label='{:4.0f} km'.format(domsize3))
        #ax.set_yscale('log')
        #plt.yticks([1000, 500, 250, 100, 50, 20])
        #ax.set_ylim(p2[-1], p2[0])
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(0, 16)
        #ax.invert_yaxis()
        #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.xlabel('{0} ({1})'.format(profile_name, profile_var1.units))
        plt.ylabel('z (km)')
        plt.title('{:s} ({:s}),  {:2.0f}-day mean at day {:3.0f} '.format(profile_name, profile_var1.units, nstat/24., t1[t1i]))
        plt.legend(loc='best')
        plt.savefig(fout + 'evoprofile_nudge_{0}days_'.format(np.round(t1[t1i]))  + profile_name + '_idealRCE.pdf')
        plt.clf()
    

#for i, ts_name in enumerate(timeseries_names):
#    ts_var1 = nc_vars1[ts_name]
#    ts1 = ts_var1[:]
#    ts_var2 = nc_vars2[ts_name]
#    ts2 = ts_var2[:]
#    ts_var3 = nc_vars3[ts_name]
#    ts3 = ts_var3[:]
#    plt.figure(1)
#    plt.plot(t1, ts1, '-k', alpha=0.8, label='{:3.0f} km'.format(domsize1))
#    plt.plot(t2, ts2, '-r', alpha=0.8, label='{:4.0f} km'.format(domsize2))
#    plt.plot(t3, ts3, '-g', alpha=0.8, label='{:4.0f} km'.format(domsize3))
#    plt.title('{0} ({1})'.format(ts_name, ts_var1.units))
#    plt.xlabel('t (days)')
#    plt.xlim((0,250))
#    plt.axvline(tdouble, color='b', alpha=0.5)
#    plt.ylabel('{0} ({1})'.format(ts_name, ts_var1.units))
#    plt.title('{:s} ({:s})'.format(ts_name, ts_var1.units))
#    plt.legend(loc='best')
#    plt.savefig(fout + 'timeseries_nudge_{0}days_'.format(np.round(t1[t3i]))  + ts_name + '_idealRCE.pdf')
#    plt.clf()
#    
#    #plot net radiation at TOA
#    if (ts_name == 'LWNTOA'):
#        plt.figure(2)
#        LWNTOA_ts1 = ts1
#        SWNTOA_ts1 = nc_vars1['SWNTOA'][:]
#        QNTOA_ts1 = SWNTOA_ts1 - LWNTOA_ts1
#        LWNTOA_ts2 = ts2
#        SWNTOA_ts2 = nc_vars2['SWNTOA'][:]
#        QNTOA_ts2 = SWNTOA_ts2 - LWNTOA_ts2
#        LWNTOA_ts3 = ts3
#        SWNTOA_ts3 = nc_vars3['SWNTOA'][:]
#        QNTOA_ts3 = SWNTOA_ts3 - LWNTOA_ts3
#        plt.title('QNTOA (W/m2)')
#        plt.xlim((0, 250))
#        plt.xlabel('t (days)')
#        plt.ylabel('QNTOA (W/m2)')
#        plt.axvline(tdouble, color='b', alpha=0.5)
#        plt.plot(t1, QNTOA_ts1, 'k', alpha=0.8, label='{:3.0f} km'.format(domsize1))
#        plt.plot(t2, QNTOA_ts2, 'r', alpha=0.8, label='{:4.0f} km'.format(domsize2))
#        plt.plot(t3, QNTOA_ts3, 'g', alpha=0.8, label='{:4.0f} km'.format(domsize3))
#        plt.legend(loc='best')
#        plt.savefig(fout + 'timeseries_nudge_{0}days_'.format(np.round(t3[t3i])) + 'QNTOA' + '_idealRCE.pdf')
#        plt.clf()
#    
#    
#    
#
#
#
