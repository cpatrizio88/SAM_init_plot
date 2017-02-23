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

matplotlib.rcParams.update({'font.size': 22})
matplotlib.rcParams.update({'figure.figsize': (20, 10)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 22})

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'

nc_in1 = glob.glob(fpath + '*256x256*3000m*day230to250*302K.nc')[0]
nc_in2 = glob.glob(fpath + '*512x512*3000m*day180to195*302K.nc')[0]
nc_in3 = glob.glob(fpath + '*1024x1024*3000m*day170to180*302K.nc')[0]

#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to130_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to150_64vert_ubarzero_MOISTDRYPROFS/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_MOISTPROFS/'

c = constants()

domsizes = [768, 1536, 3072]
ncs = [nc_in1, nc_in2, nc_in3]
colors = ['k', 'r', 'g']

#domsize=768
#domsize=1536
#domsie=3072

nave=10
ntave2D=24
ntave3D=4

ts2 = np.array([-1, -1, -1])*ntave2D
ts3 = np.array([-1, -1, -1])*ntave3D

aveperiod2D = nave*ntave2D
aveperiod3D = nave*ntave3D

for k, nc_in in enumerate(ncs):
    domsize = domsizes[k]
    print 'domsize', domsize
    
    if domsize == 768:
        nave=10
    elif domsize == 1536:
        nave=10
    else:
        nave=10
        
    aveperiod2D = nave*ntave2D
    aveperiod3D = nave*ntave3D
    
    t2 = ts2[k]
    t3 = ts3[k]
    nc_data = Dataset(nc_in)
    varis3D = nc_data.variables
    x = varis3D['x'][:]
    y = varis3D['y'][:]
    z = varis3D['z'][:]
    p = varis3D['p'][:]
    p = p*1e2
    times = varis3D['time'][t3-aveperiod3D:t3]
    print 't3D', times
    daytimes=np.arange(0, times[-1]+1)
    nx = x.size
    ny = y.size
    xx, yy = np.meshgrid(x, y)
    delx=np.diff(x)[0]
    dely=np.diff(y)[0]
    W = varis3D['W'][t3-aveperiod3D:t3,:,:,:]
    T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
    Wtave = np.mean(W, axis=0)
    
    #EDIT: this critical w threshold determines where convective region is. 
    #What is a good way to choose this value? (want to neglect gravity wave perturbations)
    Wcrit=0.01
    totpoints=nx*ny
    
    rhoz = p/(c.Rd*np.mean(np.mean(T, axis=3), axis=2))
    rhoz_tave = np.mean(rhoz, axis=0)
    
    sigma = np.zeros(z.size)
    l_ms = np.zeros(z.size)
    v_ns = np.zeros(z.size-1)
    divs = np.zeros(z.size-1)
    Aup = np.zeros(z.size)
    
    db=16
    
    #block averaging
    Wtave = blockave3D(Wtave, db)
    nxblock = nx // db
    nyblock = ny // db
    totpoints = nxblock*nyblock
    
    
    for i, zlev in enumerate(z):
        Wz = Wtave[i,:,:]
        sigma[i] = len(Wz[Wz >= Wcrit])/(1.*totpoints)
        Aup[i] = sigma[i]*(domsize*1e3)**2
        l_ms[i] = np.sqrt(Aup[i]/np.pi)
        
    #print 'vertically averaged Aup (km^2)', np.mean(Aup/1e6)
        
        
    for i in np.arange(0, z.size-1):
        delz = z[i+1] - z[i]
        Wz = Wtave[i,:,:]
        Wzup = Wtave[i+1,:,:]
        Wzup_crit = Wzup[Wz > Wcrit]
        Wz_crit = Wz[Wz > Wcrit]
        #vertical mass flux in convective region
        vertmassflux = (rhoz_tave[i+1]*Wzup_crit - rhoz_tave[i]*Wz_crit)/delz
        #total vertical mass flux convergence in convective region
        divs[i] = -np.sum(vertmassflux)*delx*dely
        
        #radial wind component (normal to convective region, assuming region is ciruclar)
        v_ns[i] = divs[i]/(2*np.pi*l_ms[i]*rhoz_tave[i])
        
    sigmabar = np.mean(sigma)
    divbar = np.mean(np.abs(divs))
    l_mbar = np.mean(l_ms/1e3)
    v_nbar = np.mean(v_ns[np.isfinite(v_ns)])
        
    print 'vertically averaged sigma', sigmabar
    print 'vertically averaged l_m (km)', l_mbar
    print 'vertically averaged v_n (m/s)', v_nbar
    print 'vertically averaged mass flux (kg/s)', np.mean(divs)
    print 'vertically averaged mass flux magnitude (kg/s)', divbar
        
        
    #f, axarr = plt.subplots(3,1)
    plt.figure(1)
    ax = plt.gcf().gca()
    if db == 1:
        plt.title('convective fractional area (W$_{{c}}$ > {:2.3f} m/s)'.format(Wcrit))
    else:
        plt.title('convective fractional area (W$_{{c}}$ > {:2.3f} m/s), block-averaging over ({:2.0f} km)$^2$'.format(Wcrit, db*(np.diff(x)[0])/1e3))
    plt.plot(sigma, z/1e3, '-x', color=colors[k], label='{:d} km, averaged over day {:2.1f} to {:2.1f}'.format(domsize, times[0], times[-1]))
    ax.set_xlabel(r'$\sigma$')
    plt.legend(loc='best')
    plt.ylabel('z (km)')
    #axarr[0,].plot(sigma, z/1e3)
    #axarr[0,].set_xlabel('sigma')
    #axarr[0,].set_ylabel('z (km)')
    
    ax.set_ylim(0, 20)
    plt.savefig(fout + 'sigma_blkavg_day250_{:d}day_db{:d}.pdf'.format(nave, db))
    
    plt.figure(2)
    ax = plt.gcf().gca()
    if db == 1:
        plt.title('normal velocity along convective region edge, W$_{{c}}$ = {:4.2} m/s'.format(Wcrit))
    else:
        plt.title('normal velocity along convective region edge, W$_{{c}}$ = {:2.1f} m/s, block-averaging over ({:2.0f} km)$^2$'.format(Wcrit, db*(np.diff(x)[0])/1e3))
    plt.plot(v_ns, z[:-1]/1e3, '-x', color=colors[k], label='{:d} km, averaged over day {:2.1f} to {:2.1f}'.format(domsize, times[0], times[-1]))
    plt.xlabel('velocity (m/s)')
    plt.ylabel('z (km)')
    plt.legend(loc='best')
    ax.set_ylim(0, 20)
    #ax.set_xlim(-10, 300)
    
    plt.savefig(fout + 'vnormal_blkavg_day250_{:d}day_db{:d}.pdf'.format(nave, db))

    plt.figure(3)
    ax = plt.gcf().gca()
    if db == 1:
        plt.title(r'horizontal mass flux out of convective region, W$_{{c}}$ = {:4.2f} m/s'.format(Wcrit))
    else:
        plt.title(r'horizontal mass flux out of convective region, W$_{{c}}$ = {:4.2f} m/s, block-averaging over ({:2.0f} km)$^2$'.format(Wcrit, db*(np.diff(x)[0])/1e3))
    plt.plot(divs/divbar, z[:-1]/1e3, '-x', color=colors[k],  label='{:d} km, averaged over day {:2.1f} to {:2.1f}'.format(domsize, times[0], times[-1]))
    plt.xlabel('mass flux relative to vertically-averaged magnitude of mass flux')
    plt.ylabel('z (km)')
    plt.legend(loc='best')
    ax.set_ylim(0, 20)
    
    plt.savefig(fout + 'massflux_blkavg_day250_{:d}day_db{:d}.pdf'.format(nave, db))
    
    #vertically smooth the divergence
    
    def moving_average(a, n=3):
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n
    
    n=5
    divsmooth = moving_average(divs/divbar, n)
    
    plt.figure(4)
    ax = plt.gcf().gca()
    
    trunc = (n-1)/2
    if db == 1:
        plt.title(r'vertically smoothed horizontal mass flux out of convective region, W$_{{c}}$ = {:4.2f} m/s, window size = {:d}'.format(Wcrit, n))
    else:
        plt.title(r'vertically smoothed horizontal mass flux out of convective region, W$_{{c}}$ = {:4.2f} m/s, window size = {:d}, block-averaging over ({:2.0f} km)$^2$'.format(Wcrit, n, db*(np.diff(x)[0])/1e3))
    plt.plot(divsmooth, z[trunc:-(trunc+1)]/1e3, '-x', color=colors[k], label='{:d} km, averaged over day {:2.1f} to {:2.1f}'.format(domsize, times[0], times[-1]))
    plt.xlabel('mass flux relative to vertically-averaged magnitude of mass flux')
    plt.ylabel('z (km)')
    plt.legend(loc='best')
    ax.set_ylim(0, 20)
    
    plt.savefig(fout + 'smoothmassflux_blkavg_day250_{:d}day_db{:d}.pdf'.format(nave, db))
    
plt.figure(1)
plt.close()

plt.figure(2)
plt.close()

plt.figure(3)
plt.close()

plt.figure(4)
plt.close()
        
        
    


#axarr[1,].plot(l_ms/1e3, z/1e3)
#axarr[1,].set_xlabel('radius of convective region (km)')
#axarr[1,].set_ylabel('z (km)')

#axarr[2,].plot(Aup/(1e6), z/1e3)
#axarr[2,].set_xlabel(r'area of convective region (km$^2$)')
#axarr[2,].set_ylabel('z (km)')







    
    
    
    
    

    
    
    



