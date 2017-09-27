import matplotlib as mpl
mpl.use('Agg')
from netCDF4 import Dataset
import site
site.addsitedir('/glade/scratch/patrizio/thermolib/')
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
#import SAM_init_plot.block_fns
#import SAM_init_plot.misc_fns
#import thermolib
#from SAM_init_plot.misc_fns import radprof3D, radprof
from misc_fns import fracclusterarea
#from SAM_init_plot.block_fns import blockave2D, blockave3D, blockxysort2D
from block_fns import blockave3D
#from thermolib.constants import constants
from constants import constants
from wsat import wsat
#from thermolib.wsat import wsat

matplotlib.rcParams.update({'font.size': 28})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 2})
matplotlib.rcParams.update({'legend.fontsize': 24})
matplotlib.rcParams.update({'mathtext.fontset': 'cm'})

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath3D = '/glade/scratch/patrizio/OUT_3D_nc/'

nc_in1 = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]
nc_in2 = glob.glob(fpath3D + '*512x512*3000m*day180to195*302K.nc')[0]
nc_in3 = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]
nc_in4 = glob.glob(fpath3D + '*2048*.nc')[0]

#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to130_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to150_64vert_ubarzero_MOISTDRYPROFS/'
fout = '/Users/cpatrizio/Google Drive/MS/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_MOISTPROFS/'

fout = '/glade/scratch/patrizio/MOISTPROF_FIGS/'

c = constants()

domsizes = [768, 1536, 3072, 6144]
ncs = [nc_in1, nc_in2, nc_in3, nc_in4]
colors = ['k', 'r', 'g', 'm']

#domsize=768
#domsize=1536
#domsie=3072

nave=10
ntave2D=24
ntave3D=4

ts2 = np.array([-1, -1, -1])
ts3 = np.array([-1, -1, -1])

aveperiod2D = nave*ntave2D
aveperiod3D = nave*ntave3D

for k, nc_in in enumerate(ncs):
    domsize = domsizes[k]
    print 'domsize', domsize
    
    if domsize == 768:
        nave=8
    elif domsize == 1536:
        nave=8
    elif domsize == 3072:
        nave = 8
    else:
        nave=5
        
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
    #Wcrit=0.5
    #totpoints=nx*ny
    
    rhoz = p/(c.Rd*np.mean(np.mean(T, axis=3), axis=2))
    rhoz_tave = np.mean(rhoz, axis=0)
    
    sigma = np.zeros(z.size)
    l_ms = np.zeros(z.size)
    v_ns = np.zeros(z.size-1)
    divs = np.zeros(z.size-1)
    Aup = np.zeros(z.size)
    
    db=1
    Wcrit=0.3
    
    #block averaging
    Wtave = blockave3D(Wtave, db)
    nxblock = nx / db
    nyblock = ny / db
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
        Wzup_crit = Wzup[Wzup > Wcrit]
        Wz_crit = Wz[Wz > Wcrit]
        #vertical mass flux in convective region
        vertmassfluxup = (rhoz_tave[i+1]*Wzup_crit)
        vertmassfluxd = (rhoz_tave[i]*Wz_crit)
        #vertmassflux = (rhoz_tave[i+1]*Wzup_crit - rhoz_tave[i]*Wz_crit)/delz
        #total vertical mass flux convergence in convective region
        divs[i] = -(np.sum(vertmassfluxup)-np.sum(vertmassfluxd))*((delx*dely)/delz)
        
        #radial wind component (normal to convective region, assuming region is ciruclar)
        v_ns[i] = divs[i]/(2*np.pi*l_ms[i]*rhoz_tave[i])
        
    sigmabar = np.mean(sigma)
    divbar = np.mean(np.abs(divs))
    l_mbar = np.mean(l_ms/1e3)
    v_nbar = np.mean(v_ns[np.isfinite(v_ns)])
    
    normmassflux = divs/divbar
    
    mid = np.bitwise_and(z > 5000, z < 7000)
    mid = mid[:-1]
    normmassflux_mid = np.mean(normmassflux[mid])

    print 'vertically averaged sigma', sigmabar
    print 'vertically averaged l_m (km)', l_mbar
    print 'vertically averaged v_n (m/s)', v_nbar
    print 'vertically averaged mass flux (kg/s)', np.mean(divs)
    print 'vertically averaged mass flux magnitude (kg/s)', divbar
    print 'mnormalized mid level mass flux', normmassflux_mid
        
        
    #f, axarr = plt.subplots(3,1)
    plt.figure(1)
    ax = plt.gcf().gca()
    if db == 1:
        tt1=plt.title(r'Convective Fractional Area')
        ax.set_xlabel(r'$\sigma$', fontsize=38)
    else:
        tt1=plt.title(r'Mesoscale Convective Fractional Area')
        ax.set_xlabel(r'$\sigma_m$', fontsize=38)
    plt.plot(sigma, z/1e3, '-x', color=colors[k], label='{:d} km, averaged over day {:2.0f} to {:2.0f}'.format(domsize, times[0], times[-1]))
    tt1.set_position([0.5, 1.04])
    #ax.set_xlabel(r'$\sigma$')
    plt.legend(loc='best')
    plt.ylabel('z (km)', fontsize=34)
    #axarr[0,].plot(sigma, z/1e3)
    #axarr[0,].set_xlabel('sigma')
    #axarr[0,].set_ylabel('z (km)')
    
    ax.set_ylim(0, 20)
    plt.savefig(fout + 'sigma_blkavg_day250_{:d}day_db{:d}_Wc{:3.3f}.pdf'.format(nave, db, Wcrit))
    
    plt.figure(2)
    ax = plt.gcf().gca()
    if db == 1:
        tt1=plt.title('normal velocity along convective region edge, $w_{{c}}$ = {:2.1} m/s'.format(Wcrit))
    else:
        tt1=plt.title('normal velocity along convective region edge, $w_{{c}}$ = {:2.2f} m/s, block-averaging over ({:2.0f} km)$^2$'.format(Wcrit, db*(np.diff(x)[0])/1e3))
    tt1.set_position([0.5, 1.04])
    plt.plot(v_ns, z[:-1]/1e3, '-x', color=colors[k], label='{:d} km, averaged over day {:2.0f} to {:2.0f}'.format(domsize, times[0], times[-1]))
    plt.xlabel('velocity (m/s)')
    plt.ylabel('z (km)')
    plt.legend(loc='best')
    ax.set_ylim(0, 20)
    #ax.set_xlim(-10, 300)
    
    plt.savefig(fout + 'vnormal_blkavg_day250_{:d}day_db{:d}_Wc{:3.3f}.pdf'.format(nave, db, Wcrit))

    plt.figure(3)
    ax = plt.gcf().gca()
    plt.gcf().subplots_adjust(bottom=0.15)
    #trunc = (n-1)/2
    if db == 1:
        tt1=plt.title(r'Horizontal Mass Flux out of Convective region, $w_{{c}}$ = {:2.1f} m/s'.format(Wcrit))
    else:
        tt1=plt.title(r'Horizontal Mass Flux out of Mesoscale Convective Region, $w_{{c}}$ = {:2.2f} m/s'.format(Wcrit))
    tt1.set_position([0.5, 1.04])  
    plt.plot(divs/divbar, z[:-1]/1e3, '-x', color=colors[k],  label='{:d} km, averaged over day {:2.0f} to {:2.0f}'.format(domsize, times[0], times[-1]))
    plt.xlabel(r'$\frac{(\nabla \cdot \rho{\bf u})_{conv}}{\overline{|\nabla \cdot \rho{\bf u}|}_{conv}}$', fontsize=40)
    plt.ylabel('z (km)', fontsize=32)
    plt.tight_layout()
    if db == 16:
        plt.xlim((-10,4))
    plt.legend(loc='best')
    ax.set_ylim(0, 20)
    
    plt.savefig(fout + 'massflux_blkavg_day250_{:d}day_db{:d}_Wc{:3.3f}.pdf'.format(nave, db, Wcrit))
    
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
        tt1=plt.title(r'Horizontal Mass Flux out of Convective Region, $w_{{c}}$ = {:2.1f} m/s'.format(Wcrit))
    else:
        tt1=plt.title(r'Horizontal Mass Flux out of Mesoscale Convective Region, $w_{{c}}$ = {:2.2f} m/s'.format(Wcrit))
    tt1.set_position([0.5, 1.04])
    plt.plot(divsmooth, z[trunc:-(trunc+1)]/1e3, '-x', color=colors[k], label='{:d} km, averaged over day {:2.0f} to {:2.0f}'.format(domsize, times[0], times[-1]))
    plt.xlabel(r'$\frac{(\nabla \cdot \rho{\bf u})_{conv}}{\overline{|\nabla \cdot \rho{\bf u}|}_{conv}}$', fontsize=40)
    plt.ylabel('z (km)')
    plt.tight_layout()
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







    
    
    
    
    

    
    
    



