import matplotlib as mpl
mpl.use('Agg')
from netCDF4 import Dataset
import site
site.addsitedir('/glade/scratch/patrizio/thermolib/')
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker
import matplotlib
import numpy.ma as ma
#from thermolib.constants import constants
from constants import constants
from wsat import wsat
#from thermolib.wsat import wsat

c = constants()

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

matplotlib.rcParams.update({'font.size': 28})
#matplotlib.rcParams.update({'figure.figsize': (18, 10)})
matplotlib.rcParams.update({'lines.linewidth': 1})
matplotlib.rcParams.update({'legend.fontsize': 22})

matplotlib.rcParams.update({'mathtext.fontset': 'cm'})

plt.style.use('seaborn-white')

fpath = '/glade/scratch/patrizio/OUT_STAT_nc/'
#fout = '/Users/cpatrizio/Google Drive/MS/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STAT/'
fout = '/glade/scratch/patrizio/OUT_STAT_FIGS/'

nc_in1 = glob.glob(fpath + '*256x256*3000m*250days_302K.nc')[0]
nc_in2 = glob.glob(fpath + '*512x512*3000m*195days_302K.nc')[0]
nc_in3 = glob.glob(fpath + '*1024x1024*3000m*220days_302K.nc')[0]
nc_in4 = glob.glob(fpath + '*2048x2048*3000m*300days_302K.nc')[0]

delx1=3000
domsize1 = 768
delx2=3000
domsize2 = 1536
domsize3 = 3072
domsize4 = 6144

tdouble=90

nc_data1 = Dataset(nc_in1)
nc_data2 = Dataset(nc_in2)
nc_data3 = Dataset(nc_in3)
nc_data4 = Dataset(nc_in4)
nc_vars1 = nc_data1.variables
nc_vars2 = nc_data2.variables
nc_vars3 = nc_data3.variables
nc_vars4 = nc_data4.variables



z1= nc_vars1['z'][:]
z2 = nc_vars2['z'][:]
z3 = nc_vars3['z'][:]
z4 = nc_vars4['z'][:]

p1 = nc_vars1['p'][:]
p2 = nc_vars2['p'][:]
p3 = nc_vars3['p'][:]
p4 = nc_vars4['p'][:]

t1 = nc_vars1['time'][:]
t2 = nc_vars2['time'][:]
t3 = nc_vars3['time'][:]


t4d = -24*97 #time when doubling occurs

t4 = nc_vars4['time'][t4d:]


tt1, zz1 = np.meshgrid(t1, z1)
tt1, pp1 = np.meshgrid(t1, p1)

dt1 = (t1[-1] - t1[-2]) #difference between OUT_STAT output in days
dt2 = (t2[-1] - t2[-2]) 
dt3 = (t3[-1] - t3[-2])


tt2, zz2 = np.meshgrid(t2, z2)
tt2, pp2 = np.meshgrid(t2, p2)

tt3, zz3 = np.meshgrid(t3, z3)
tt3, pp3 = np.meshgrid(t3, p3)

tt4, zz4 = np.meshgrid(t4, z4)
tt4, pp4 = np.meshgrid(t4, p4)

profile_names = ['RELH', 'PRECIP', 'THETAE', 'THETA', 'QN', 'QV', 'QCL', 'QP', 'QPI', 'QCOND', 'TABS', 'U', 'V', 'MSE', 'DSE', 'RADQR', ]

#profile_names = ['QI', 'CAPE']
timeseries_names = ['PW', 'SHF', 'PREC', 'SWNTOA', 'SWNS', 'CLDLOW', 'CLDHI', 'LWNTOA','LWNS', 'LHF',] 

profile_names = ['QPI']

#profile_names = ['QP']

#profile_names = ['QPL']

nstat=24

taggr = -1
t2i = -1
t3i = -1
t4i = -1

nave = int(10*nstat) #average the profiles over nave hours to get the daily time tendency

for i, profile_name in enumerate(profile_names):
    
    if profile_name == 'PRECIP':
        titlename='P'
    elif profile_name == 'RELH':
        titlename='RH'
    elif profile_name == 'THETAE':
        titlename=r'$\theta_e$'
    elif profile_name == 'THETA':
        titlename=r'$\theta$'
    elif profile_name == 'QN':
        titlename=r'$q_n$'
    elif profile_name == 'QI':
        titlename = r'$q_i$'
    elif profile_name == 'QV':
        titlename = r'$q_v$'
    elif profile_name == 'QCL':
        titlename = r'$q_{cl}$'
    elif profile_name == 'QP':
        titlename = r'$q_p$'
    elif profile_name == 'QPI':
        titlename = r'$q_{pi}$'
    elif profile_name == 'QCOND':
        titlename = r'$q_p + q_n$'
    elif profile_name == 'TABS':
        titlename = 'T'
    elif profile_name == 'RADQR':
        titlename = r'$Q_r$'
    else:
        titlename = profile_name
    
    profile_var1 = nc_vars1[profile_name]
    profile1 = profile_var1[:]

    profile_var2 = nc_vars2[profile_name]
    profile2 = profile_var2[:]
    
    profile_var3 = nc_vars3[profile_name]
    profile3 = profile_var3[:]
    
    profile_var4 = nc_vars4[profile_name]
    profile4 = profile_var4[t4d:,:]
    
    tendaveprofile1 = np.mean(profile1[taggr-nave:taggr,:], axis=0)
    #tendaveprofile_prev1 = np.sum(profile1[-2*nave:-nave, :], axis=0)/nave
 
    tendaveprofile2 = np.mean(profile2[t2i-nave:t2i,:], axis=0)
    #tendaveprofile_prev2 = np.sum(profile2[-2*nave:-nave, :], axis=0)/nave
    
    tendaveprofile3 = np.mean(profile3[t3i-nave:t3i,:], axis=0)
    
    tendaveprofile4 = np.mean(profile4[t4i-nave:t4i,:], axis=0)
    
    if profile_name == 'QN':
        zs1 = np.zeros(profile1.shape)
        zs2 = np.zeros(profile2.shape)
        zs3 = np.zeros(profile3.shape)
        zs4 = np.zeros(profile4.shape)
        zs1[:,:] = z1[::-1]
        zs2[:,:] = z2[::-1]
        zs3[:,:] = z3[::-1]
        zs4[:,:] = z4[::-1]
        cloudthresh = 0.001 #threshold cumulative condensate 
        profile_cumsum1 = np.cumsum(profile1[:,::-1], axis=1)
        #cloudtop1 = np.max(zs[profile_cumsum1 > cloudthresh], axis=1)
        profile_cumsum2 = np.cumsum(profile2[:,::-1], axis=1)
        #cloudtop2 = np.max(zs[profile_cumsum2 > cloudthresh], axis=1)
        profile_cumsum3 = np.cumsum(profile3[:,::-1], axis=1)
        #cloudtop3 = np.max(zs[profile_cumsum3 > cloudthresh], axis=1)
        profile_cumsum4 = np.cumsum(profile4[:,::-1], axis=1)
    
        cloudtop1 = ma.masked_array(zs1, mask=profile_cumsum1 < cloudthresh)
        cloudtop1 = np.max(cloudtop1, axis=1)
        cloudtop2 = ma.masked_array(zs2, mask=profile_cumsum2 < cloudthresh)
        cloudtop2 = np.max(cloudtop2, axis=1)
        cloudtop3 = ma.masked_array(zs3, mask=profile_cumsum3 < cloudthresh)
        cloudtop3 = np.max(cloudtop3, axis=1)
        cloudtop4 = ma.masked_array(zs4, mask=profile_cumsum4 < cloudthresh)
        cloudtop4 = np.max(cloudtop4, axis=1)
        
        
    if profile_name == 'QPI':
        T1 = nc_vars1['TABS'][:]
        p2D1 = np.zeros(T1.shape)
        p2D1[:,:] = p1*1e2
        rho1 = (c.Rd*T1)/p2D1
        T2 = nc_vars2['TABS'][:]
        p2D2 = np.zeros(T2.shape)
        p2D2[:,:] = p2*1e2
        rho2 = (c.Rd*T2)/p2D2
        T3 = nc_vars3['TABS'][:]
        p2D3 = np.zeros(T3.shape)
        p2D3[:,:] = p3*1e2
        rho3 = (c.Rd*T3)/p2D3
        T4 = nc_vars3['TABS'][t4d:,:]
        p2D4 = np.zeros(T4.shape)
        p2D4[:,:] = p4*1e2
        rho4 = (c.Rd*T4)/p2D4
  
        
        diffz1 = np.diff(z1)
        diffz2D1 = np.zeros((T1.shape[0], T1.shape[1]-1))
        diffz2D1[:,:] = diffz1
        diffz2 = np.diff(z2)
        diffz2D2 = np.zeros((T2.shape[0], T2.shape[1]-1))
        diffz2D2[:,:] = diffz2
        diffz3 = np.diff(z3)
        diffz2D3 = np.zeros((T3.shape[0], T3.shape[1]-1))
        diffz2D3[:,:] = diffz3
        diffz4 = np.diff(z4)
        diffz2D4 = np.zeros((T4.shape[0], T4.shape[1]-1))
        diffz2D4[:,:] = diffz4

        PIpath1 = np.sum(rho1[:,:-1]*profile1[:,:-1]*diffz2D1, axis=1)
        PIpath2 = np.sum(rho2[:,:-1]*profile2[:,:-1]*diffz2D2, axis=1)
        PIpath3 = np.sum(rho3[:,:-1]*profile3[:,:-1]*diffz2D3, axis=1)
        PIpath4 = np.sum(rho4[:,:-1]*profile4[:,:-1]*diffz2D4, axis=1)
#        
#    if profile_name == 'QP':
#        T1 = nc_vars1['TABS'][:]
#        p2D1 = np.zeros(T1.shape)
#        p2D1[:,:] = p1*1e2
#        rho1 = (c.Rd*T1)/p2D1
#        T2 = nc_vars2['TABS'][:]
#        p2D2 = np.zeros(T2.shape)
#        p2D2[:,:] = p2*1e2
#        rho2 = (c.Rd*T2)/p2D2
#        T3 = nc_vars3['TABS'][:]
#        p2D3 = np.zeros(T3.shape)
#        p2D3[:,:] = p3*1e2
#        rho3 = (c.Rd*T3)/p2D3
#  
#        
#        diffz1 = np.diff(z1)
#        diffz2D1 = np.zeros((T1.shape[0], T1.shape[1]-1))
#        diffz2D1[:,:] = diffz1
#        diffz2 = np.diff(z2)
#        diffz2D2 = np.zeros((T2.shape[0], T2.shape[1]-1))
#        diffz2D2[:,:] = diffz2
#        diffz3 = np.diff(z3)
#        diffz2D3 = np.zeros((T3.shape[0], T3.shape[1]-1))
#        diffz2D3[:,:] = diffz3
#
#        PIpath1 = np.sum(rho1[:,:-1]*profile1[:,:-1]*diffz2D1, axis=1)
#        PIpath2 = np.sum(rho2[:,:-1]*profile2[:,:-1]*diffz2D2, axis=1)
#        PIpath3 = np.sum(rho3[:,:-1]*profile3[:,:-1]*diffz2D3, axis=1)
#        
#    if profile_name == 'QPL':
#        T1 = nc_vars1['TABS'][:]
#        p2D1 = np.zeros(T1.shape)
#        p2D1[:,:] = p1*1e2
#        rho1 = (c.Rd*T1)/p2D1
#        T2 = nc_vars2['TABS'][:]
#        p2D2 = np.zeros(T2.shape)
#        p2D2[:,:] = p2*1e2
#        rho2 = (c.Rd*T2)/p2D2
#        T3 = nc_vars3['TABS'][:]
#        p2D3 = np.zeros(T3.shape)
#        p2D3[:,:] = p3*1e2
#        rho3 = (c.Rd*T3)/p2D3
#  
#        
#        diffz1 = np.diff(z1)
#        diffz2D1 = np.zeros((T1.shape[0], T1.shape[1]-1))
#        diffz2D1[:,:] = diffz1
#        diffz2 = np.diff(z2)
#        diffz2D2 = np.zeros((T2.shape[0], T2.shape[1]-1))
#        diffz2D2[:,:] = diffz2
#        diffz3 = np.diff(z3)
#        diffz2D3 = np.zeros((T3.shape[0], T3.shape[1]-1))
#        diffz2D3[:,:] = diffz3
#
#        PIpath1 = np.sum(rho1[:,:-1]*profile1[:,:-1]*diffz2D1, axis=1)
#        PIpath2 = np.sum(rho2[:,:-1]*profile2[:,:-1]*diffz2D2, axis=1)
#        PIpath3 = np.sum(rho3[:,:-1]*profile3[:,:-1]*diffz2D3, axis=1)
        
        
        
        
        
        
 
    #difference between model start and model end profiles
    plt.figure(1, figsize=(18,12))
    ax = plt.gca()
    index_mid = np.where(t1 >= 60)[0][0]
    #t0profile1 = profile1[0,:]
    #tmidprofile1 = profile1[index_mid, :]
    #tendprofile1 = profile1[-1,:]
    #t0profile2 = profile2[0,:]
    #tmidprofile2 = profile2[index_mid, :]
    #tendprofile2 = profile2[-1,:]
    #plt.plot(t0profile1, p1, 'k--', label='t = 0')
    #plt.plot(tmidprofile1, p1, 'b-', label='{0} at t = {1} days'.format(profile_name, np.round(t1[index_mid])))
    #plt.plot(tendaveprofile_prev1, p1, 'r-', label='{:2.0f}-day time average t = {:3.1f} days'.format(nave/24., t1[-nave]))
    plt.plot(tendaveprofile1, z1/1e3, 'k-', alpha=0.8, linewidth=2.5, label='{:3.0f} km, {:2.0f}-day time average at day {:3.0f}'.format(domsize1, nave/24., t1[taggr]))
    #plt.plot(t0profile2, p2, 'k--', alpha=0.2, label='ubar = 0, {0} at t = 0'.format(profile_name))
    #plt.plot(tmidprofile2, p2, 'b-', alpha=0.2, label='ubar = 0, {0} at t = {1} days'.format(profile_name, np.round(t1[index_mid])))
    #plt.plot(tendaveprofile_prev2, p2, 'r-', alpha=0.2, label='{:2.0f}-day time average at t = {:3.1f} days'.format(nave/24.,t1[-nave]))
    plt.plot(tendaveprofile2, z2/1e3, 'r-', alpha=0.8, linewidth=2.5, label='{:4.0f} km, {:2.0f}-day time average at day {:3.0f}'.format(domsize2, nave/24., t2[t2i]))
    plt.plot(tendaveprofile3, z3/1e3, 'g-', alpha=0.8, linewidth=2.5, label='{:4.0f} km, {:2.0f}-day time average at day {:3.0f}'.format(domsize3, nave/24., t3[t3i]))
    plt.plot(tendaveprofile4, z4/1e3, 'm-', alpha=0.8, linewidth=2.5, label='{:4.0f} km, {:2.0f}-day time average at day {:3.0f}'.format(domsize4, nave/24., t4[t4i]))
    #ax.set_yscale('log')
    #plt.yticks([1000, 500, 250, 100, 50, 20])
    #ax.set_ylim(p2[-1], p2[0])
    #ax.invert_yaxis()
    #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xlabel('{0} ({1})'.format(titlename, profile_var1.units))
    plt.ylabel('z (km)')
    plt.ylim(0, 20)
    plt.title('{:s}'.format(titlename))
    plt.legend(loc='best')
    plt.savefig(fout + 'evoprofile_nudge_{0}days_'.format(np.round(t1[t3i]))  + profile_name + '_idealRCE.pdf')
    plt.close()
    
    
    

for i, ts_name in enumerate(timeseries_names):

#timeseries_names = ['PW', 'LHF', 'SHF', 'PREC', 'CAPE', 'LWNTOA', 'LWNT', 'SWNTOA', 'LWDS', 'SWNS', 'CLDLOW', 'CLDHI'] 

    if ts_name == 'PREC':
        titlename='P'
    elif ts_name == 'LWNTOA':
        titlename='Q_{net,TOA}'
    elif ts_name == 'SWNTOA':
        titlename='Q_{net,TOA}'
    elif ts_name == 'SWNS':
        titlename='$Q_{net,s}'
    elif ts_name == 'LWNS':
        titlename='Q_{net,s}'
    elif ts_name == 'CLDLOW':
        titlename = 'f_{c,low}'
    elif ts_name == 'CLDHI':
        titlename = 'f_{c,high}'
    else:
        titlename = ts_name
    
    
    plt.figure(1, figsize=(18,10))
    ts_var1 = nc_vars1[ts_name]
    ts1 = ts_var1[:]
    ts_var2 = nc_vars2[ts_name]
    ts2 = ts_var2[:]
    ts_var3 = nc_vars3[ts_name]
    ts3 = ts_var3[:]
    ts_var4 = nc_vars4[ts_name]
    ts4 = ts_var4[:]
    
    if ts_name == 'PREC' or 'LHF':
        w = nstat*5
        ts4_smooth = moving_average(ts4[t4d:], w)
        ts3_smooth = moving_average(ts3, w)
        ts2_smooth = moving_average(ts2, w)
        ts1_smooth = moving_average(ts1, w)
        plt.plot(t1, ts1, '-k', alpha=0.5, label='{:3.0f} km'.format(domsize1))
        plt.plot(t2, ts2, '-r', alpha=0.5, label='{:4.0f} km'.format(domsize2))
        plt.plot(t3, ts3, '-g', alpha=0.5, label='{:4.0f} km'.format(domsize3))
        plt.plot(t4[t4d:], ts4[t4d:], '-m', alpha=0.5, label='{:4.0f} km'.format(domsize4))
        plt.plot(t1[w/2:-w/2+1], ts1_smooth, '-k', linewidth=2)
        plt.plot(t2[w/2:-w/2+1], ts2_smooth, '-r', linewidth=2)
        plt.plot(t3[w/2:-w/2+1], ts3_smooth, '-g', linewidth=2)
        plt.plot(t4[t4d+w/2:-w/2+1], ts4_smooth, '-m', linewidth=2)
    else:
        plt.plot(t1, ts1, '-k', label='{:3.0f} km'.format(domsize1))
        plt.plot(t2, ts2, '-r', label='{:4.0f} km'.format(domsize2))
        plt.plot(t3, ts3, '-g', label='{:4.0f} km'.format(domsize3))
        plt.plot(t4[t4d:], ts4[t4d:], '-m', label='{:3.0f} km'.format(domsize4))

    #plt.title('{0} ({1})'.format(ts_name, ts_var1.units))
    plt.xlabel('time (days)')
    plt.xlim((0,300))
    plt.axvline(tdouble, color='b', alpha=0.5)
    if titlename == 'PW' or titlename == 'LHF' or titlename == 'SHF' or titlename == 'P':
        plt.ylabel(r'$\overline{{\mathrm{{{0}}}}}$ ({1})'.format(titlename, ts_var1.units))
        plt.title(r'$\overline{{\mathrm{{{:s}}}}}$'.format(titlename))
    else:
        plt.ylabel(r'$\overline{{{0}}}$ ({1})'.format(titlename, ts_var1.units))
        plt.title(r'$\overline{{{0}}}$'.format(titlename)) 
    #else:
    #plt.title('{:s} ({:s})'.format(ts_name, ts_var1.units))
    plt.legend(loc='best')
    plt.savefig(fout + 'LARGEtimeseries_nudge_{0}days_'.format(np.round(t1[t3i]))  + ts_name + '_idealRCE.pdf')
    plt.close()
    
    #plot net radiation at TOA
    if (ts_name == 'LWNTOA'):
        plt.figure(2, figsize=(18,10))
        LWNTOA_ts1 = ts1
        SWNTOA_ts1 = nc_vars1['SWNTOA'][:]
        QNTOA_ts1 = SWNTOA_ts1 - LWNTOA_ts1
        LWNTOA_ts2 = ts2
        SWNTOA_ts2 = nc_vars2['SWNTOA'][:]
        QNTOA_ts2 = SWNTOA_ts2 - LWNTOA_ts2
        LWNTOA_ts3 = ts3
        SWNTOA_ts3 = nc_vars3['SWNTOA'][:]
        QNTOA_ts3 = SWNTOA_ts3 - LWNTOA_ts3
        LWNTOA_ts4 = ts4
        SWNTOA_ts4 = nc_vars4['SWNTOA'][:]
        QNTOA_ts4 = SWNTOA_ts4 - LWNTOA_ts4
        plt.title(r'$Q_{net,TOA}$')
        plt.xlim((0, 276))
        plt.xlabel('time (days)')
        plt.ylabel(r'$Q_{net,TOA}$ (W/m$^2$)')
        plt.axvline(tdouble, color='b', alpha=0.5)
        plt.plot(t1, QNTOA_ts1, 'k', alpha=0.8, label='{:3.0f} km'.format(domsize1))
        plt.plot(t2, QNTOA_ts2, 'r', alpha=0.8, label='{:4.0f} km'.format(domsize2))
        plt.plot(t3, QNTOA_ts3, 'g', alpha=0.8, label='{:4.0f} km'.format(domsize3))
        plt.plot(t4[t4d:], QNTOA_ts4[t4d:], 'm', alpha=0.8, label='{:4.0f} km'.format(domsize4))
        plt.legend(loc='best')
        plt.savefig(fout + 'LARGEtimeseries_nudge_{0}days_'.format(np.round(t1[t3i])) + 'QNTOA' + '_idealRCE.pdf')
        plt.close()
        
        plt.figure()
        LWN_ts1 = nc_vars1['LWNS'][:] - LWNTOA_ts1
        LWN_ts2 = nc_vars2['LWNS'][:] - LWNTOA_ts2
        LWN_ts3 = nc_vars3['LWNS'][:] - LWNTOA_ts3
        LWN_ts4 = nc_vars4['LWNS'][:] - LWNTOA_ts4
        plt.title(r'$Q_{LW,net}$')
        plt.xlim((0, 300))
        plt.xlabel('time (days)')
        plt.ylabel(r'$Q_{LW,net}$ (W/m$^2$)')
        plt.axvline(tdouble, color='b', alpha=0.5)
        plt.plot(t1, LWN_ts1, 'k', alpha=0.8, label='{:3.0f} km'.format(domsize1))
        plt.plot(t2, LWN_ts2, 'r', alpha=0.8, label='{:4.0f} km'.format(domsize2))
        plt.plot(t3, LWN_ts3, 'g', alpha=0.8, label='{:4.0f} km'.format(domsize3))
        plt.plot(t4[t4d:], LWN_ts4[t4d:], 'm', alpha=0.8, label='{:4.0f} km'.format(domsize4))
        plt.legend(loc='best')
        plt.savefig(fout + 'LARGEtimeseries_nudge_{0}days_'.format(np.round(t1[t3i])) + 'LWN' + '_idealRCE.pdf')
        plt.close()
        
        plt.figure()
        SWN_ts1 = SWNTOA_ts1 - nc_vars1['SWNS'][:] 
        SWN_ts2 = SWNTOA_ts2 - nc_vars2['SWNS'][:] 
        SWN_ts3 = SWNTOA_ts3 - nc_vars3['SWNS'][:] 
        SWN_ts4 = SWNTOA_ts4 - nc_vars4['SWNS'][:] 
        plt.title(r'$Q_{SW,net}$')
        plt.xlim((0, 300))
        plt.xlabel('time (days)')
        plt.ylabel(r'$Q_{SW,net}$ (W/m$^2$)')
        plt.axvline(tdouble, color='b', alpha=0.5)
        plt.plot(t1, SWN_ts1, 'k', alpha=0.8, label='{:3.0f} km'.format(domsize1))
        plt.plot(t2, SWN_ts2, 'r', alpha=0.8, label='{:4.0f} km'.format(domsize2))
        plt.plot(t3, SWN_ts3, 'g', alpha=0.8, label='{:4.0f} km'.format(domsize3))
        plt.plot(t4[t4d:], SWN_ts4[t4d:], 'm', alpha=0.8, label='{:4.0f} km'.format(domsize4))
        plt.legend(loc='best')
        plt.savefig(fout + 'LARGEtimeseries_nudge_{0}days_'.format(np.round(t1[t3i])) + 'SWN' + '_idealRCE.pdf')
        plt.close()
        
   
    if ts_name == 'LWNS':        
        plt.figure(3, figsize=(18,10))
        LWNS_ts1 = ts1
        SWNS_ts1 = nc_vars1['SWNS'][:]
        QNS_ts1 = SWNS_ts1 - LWNS_ts1
        LWNS_ts2 = ts2
        SWNS_ts2 = nc_vars2['SWNS'][:]
        QNS_ts2 = SWNS_ts2 - LWNS_ts2
        LWNS_ts3 = ts3
        SWNS_ts3 = nc_vars3['SWNS'][:]
        QNS_ts3 = SWNS_ts3 - LWNS_ts3
        LWNS_ts4 = ts4
        SWNS_ts4 = nc_vars4['SWNS'][:]
        QNS_ts4 = SWNS_ts4 - LWNS_ts4
        QNATMOS_ts1 = QNTOA_ts1 - QNS_ts1
        QNATMOS_ts2 = QNTOA_ts2 - QNS_ts2
        QNATMOS_ts3 = QNTOA_ts3 - QNS_ts3
        QNATMOS_ts4 = QNTOA_ts4 - QNS_ts4
        plt.title(r'$Q_{net,atmos}$')
        plt.xlim((0, 300))
        plt.xlabel('time (days)')
        plt.ylabel(r'$Q_{net,atmos}$ (W/m$^2$)')
        plt.axvline(tdouble, color='b', alpha=0.5)
        plt.plot(t1, QNATMOS_ts1, 'k', alpha=0.8, label='{:3.0f} km'.format(domsize1))
        plt.plot(t2, QNATMOS_ts2, 'r', alpha=0.8, label='{:4.0f} km'.format(domsize2))
        plt.plot(t3, QNATMOS_ts3, 'g', alpha=0.8, label='{:4.0f} km'.format(domsize3))
        plt.plot(t4[t4d:], QNATMOS_ts4[t4d:], 'm', alpha=0.8, label='{:4.0f} km'.format(domsize4))
        plt.legend(loc='best')
        plt.savefig(fout + 'LARGEtimeseries_nudge_{0}days_'.format(np.round(t1[t3i])) + 'QNATMOS' + '_idealRCE.pdf')
        plt.close()
        
    if ts_name == 'LHF': 
        plt.figure(4, figsize=(18,10))
        SHF_ts1 = nc_vars1['SHF'][:]
        SHF_ts2 = nc_vars2['SHF'][:]
        SHF_ts3 = nc_vars3['SHF'][:]
        SHF_ts4 = nc_vars3['SHF'][:]
        NRGBALANCE_ts1 = ts1 + QNATMOS_ts1
        NRGBALANCE_ts2 = ts2 + QNATMOS_ts2
        NRGBALANCE_ts3 = ts3 + QNATMOS_ts3 
        NRGBALANCE_ts4 = ts4 + QNATMOS_ts4 
        plt.title(r'LHF + $Q_{net,atmos}$')
        plt.xlim((0, 276))
        plt.xlabel('time (days)')
        plt.ylabel(r'LHF  + $Q_{net,atmos}$ (W/m$^2$)')
        plt.axvline(tdouble, color='b', alpha=0.5)
        plt.plot(t1, NRGBALANCE_ts1, 'k', alpha=0.8, label='{:3.0f} km'.format(domsize1))
        plt.plot(t2, NRGBALANCE_ts2, 'r', alpha=0.8, label='{:4.0f} km'.format(domsize2))
        plt.plot(t3, NRGBALANCE_ts3, 'g', alpha=0.8, label='{:4.0f} km'.format(domsize3))
        plt.plot(t4[t4d:], NRGBALANCE_ts4[t4d:], 'm', alpha=0.8, label='{:4.0f} km'.format(domsize4))
        plt.legend(loc='best')
        plt.savefig(fout + 'LARGEtimeseries_nudge_{0}days_'.format(np.round(t1[t3i])) + 'LHFQbalance' + '_idealRCE.pdf')
        plt.close()

              
#plt.figure(5, figsize=(18,10))
#plt.axvline(tdouble, color='b', alpha=0.5)
#plt.plot(t1, cloudtop1/1e3, 'k', alpha=0.8, label='{:3.0f} km'.format(domsize1))
#plt.plot(t2, cloudtop2/1e3, 'r', alpha=0.8, label='{:4.0f} km'.format(domsize2))
#plt.plot(t3, cloudtop3/1e3, 'g', alpha=0.8, label='{:4.0f} km'.format(domsize3))
#plt.title('max cloud top height')
#plt.xlabel('time (days)')
#plt.ylabel('max cloud top height (km)')
#plt.legend(loc='best')
#plt.savefig(fout + 'timeseries_nudge_{0}days_'.format(np.round(t1[t3i])) + 'cloudtop' + '_idealRCE.pdf')
#plt.clf()


PIpath1_smooth = moving_average(PIpath1/1e3, w)
PIpath2_smooth = moving_average(PIpath2/1e3, w)
PIpath3_smooth = moving_average(PIpath3/1e3, w)
PIpath4_smooth = moving_average(PIpath4/1e3, w)

plt.figure(6, figsize=(18,10))
plt.clf()
plt.axvline(tdouble, color='b', alpha=0.8)
plt.plot(t1, PIpath1/1e3, 'k', alpha=0.5, label='{:3.0f} km'.format(domsize1))
plt.plot(t2[tdouble*nstat:], PIpath2[tdouble*nstat:]/1e3, 'r', alpha=0.5, label='{:4.0f} km'.format(domsize2))
plt.plot(t3[tdouble*nstat:], PIpath3[tdouble*nstat:]/1e3, 'g', alpha=0.5, label='{:4.0f} km'.format(domsize3))
plt.plot(t4, PIpath4/1e3, 'm', alpha=0.5, label='{:4.0f} km'.format(domsize4))
plt.plot(t1[w/2:-w/2+1], PIpath1_smooth, '-k', linewidth=2.5)
plt.plot(t2[tdouble*nstat:-w/2+1], PIpath2_smooth[tdouble*nstat-w/2:], '-r', linewidth=2.5)
plt.plot(t3[tdouble*nstat:-w/2+1], PIpath3_smooth[tdouble*nstat-w/2:], '-g', linewidth=2.5)
plt.plot(t4[w/2:-w/2+1], PIpath4_smooth, 'm', linewidth=2.5)
plt.title(r'Precipitating Ice Path, $\int \rho q_{p,i} dz$')
plt.xlabel('time (days)')
plt.ylim((0, 0.2))
plt.ylabel(r'$\int \rho q_{p,i} dz$ (kg m$^{-2}$)')
plt.legend(loc='best')
plt.savefig(fout + 'timeseries_nudge_{0}days_'.format(np.round(t1[t3i])) + 'intQPI' + '_idealRCE.pdf')
plt.close()

nave=40

PIpath1bar = np.mean(PIpath1[-nave*24:])/1e3
PIpath2bar = np.mean(PIpath2[-nave*24:])/1e3
PIpath3bar = np.mean(PIpath3[-nave*24:])/1e3
PIpath4bar = np.mean(PIpath4[-nave*24:])/1e3

print 'average q_p,i, 768 km domain', PIpath1bar
print 'average q_p,i, 1536 km domain', PIpath2bar
print 'average q_p,i, 3072 km domain', PIpath3bar
print 'average q_p,i, 6144 km domain', PIpath4bar

# plt.figure(6, figsize=(18,10))
# plt.clf()
# plt.axvline(tdouble, color='b', alpha=0.5)
# plt.plot(t1, PIpath1/1e3, 'k', alpha=0.8, label='{:3.0f} km'.format(domsize1))
# plt.plot(t2, PIpath2/1e3, 'r', alpha=0.8, label='{:4.0f} km'.format(domsize2))
# plt.plot(t3, PIpath3/1e3, 'g', alpha=0.8, label='{:4.0f} km'.format(domsize3))
# plt.plot(t4, PIpath4/1e3, 'm', alpha=0.8, label='{:4.0f} km'.format(domsize4))
# plt.title(r'Ice Precipitation Path, $\int \rho Q_{p,i} dz$')
# plt.xlabel('time (days)')
# plt.ylim(0, 0.3)
# plt.ylabel(r'$\int \rho q_{p,i} dz$ (kg m$^{-2}$)')
# plt.legend(loc='best')
# plt.savefig(fout + 'timeseries_nudge_{0}days_'.format(np.round(t1[t3i])) + 'intQPI' + '_idealRCE.pdf')
# plt.close()

#plt.figure(6, figsize=(18,10))
#plt.clf()
#plt.axvline(tdouble, color='b', alpha=0.5)
#plt.plot(t1, PIpath1/1e3, 'k', alpha=0.8, label='{:3.0f} km'.format(domsize1))
#plt.plot(t2, PIpath2/1e3, 'r', alpha=0.8, label='{:4.0f} km'.format(domsize2))
#plt.plot(t3, PIpath3/1e3, 'g', alpha=0.8, label='{:4.0f} km'.format(domsize3))
#plt.plot(t4, PIpath4/1e3, 'm', alpha=0.8, label='{:4.0f} km'.format(domsize4))
#plt.title(r'liquid precipitation path, $\int \rho Q_{p,l} dz$')
#plt.xlabel('time (days)')
#plt.ylabel(r'$\int \rho Q_{p,l} dz$ (kg m$^{-2}$)')
#plt.legend(loc='best')
#plt.savefig(fout + 'timeseries_nudge_{0}days_'.format(np.round(t1[t3i])) + 'intQPL' + '_idealRCE.pdf')
#plt.close()

#plt.figure(6, figsize=(18,10))
#plt.clf()
#plt.axvline(tdouble, color='b', alpha=0.5)
#plt.plot(t1, PIpath1/1e3, 'k', alpha=0.8, label='{:3.0f} km'.format(domsize1))
#plt.plot(t2, PIpath2/1e3, 'r', alpha=0.8, label='{:4.0f} km'.format(domsize2))
#plt.plot(t3, PIpath3/1e3, 'g', alpha=0.8, label='{:4.0f} km'.format(domsize3))
#plt.title(r'total precipitation path, $\int \rho Q_{p} dz$')
#plt.xlabel('time (days)')
#plt.ylabel(r'$\int \rho Q_{p} dz$ (kg m$^{-2}$)')
#plt.legend(loc='best')
#plt.savefig(fout + 'timeseries_nudge_{0}days_'.format(np.round(t1[t3i])) + 'intQP' + '_idealRCE.pdf')
#plt.close()
        
        
    
    
    



