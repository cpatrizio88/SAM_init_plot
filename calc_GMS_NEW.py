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

matplotlib.rcParams.update({'font.size': 26})
matplotlib.rcParams.update({'figure.figsize': (20, 12)})
matplotlib.rcParams.update({'lines.linewidth': 2})
matplotlib.rcParams.update({'legend.fontsize': 24})
plt.style.use('seaborn-white')

g=9.81

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'

nc_in1 = glob.glob(fpath + '*256x256*3000m*day230to250*302K.nc')[0]
nc_in2 = glob.glob(fpath + '*512x512*3000m*day180to195*302K.nc')[0]
nc_in3 = glob.glob(fpath + '*1024x1024*3000m*day170to180*302K.nc')[0]

#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to130_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to150_64vert_ubarzero_MOISTDRYPROFS/'
fout = '/Users/cpatrizio/Google Drive/MS/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_MOISTPROFS/'

c = constants()
L_v = 2.257*1e6

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
        nave=5
    elif domsize == 1536:
        nave=5
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
    dt = np.diff(times)[0]
    dz = np.mean(np.diff(z))
    print 't3D', times
    #daytimes=np.arange(0, times[-1]+1)
    nx = x.size
    ny = y.size
    nz = z.size
    nt = times.size
    xx, yy = np.meshgrid(x, y)
    delx=np.diff(x)[0]
    dely=np.diff(y)[0]
    W = varis3D['W'][t3-aveperiod3D:t3,:,:,:]
    U = varis3D['U'][t3-aveperiod3D:t3,:,:,:]
    V = varis3D['V'][t3-aveperiod3D:t3,:,:,:]
    T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
    QV = varis3D['QV'][t3-aveperiod3D:t3,:,:,:]
    QV = QV*1e-3
    Wtave = np.mean(W, axis=0)
    Utave = np.mean(U, axis=0)
    Vtave = np.mean(V, axis=0)
    z3D = np.zeros((ny, nx, nz))
    z3D[:,:,:] = z
    z3D = z3D[:,:,:,np.newaxis]
    z4D = np.tile(z3D,nt)
    z4D = z4D.T
    z3D = np.zeros(Utave.shape).T
    #z3D[:,:,:] = z
    #z3D = z3D.T
    s = c.cpd*T + g*z4D
    h = s + L_v*QV
    hu = h*U
    hv = h*V
    
    dhudx = np.gradient(hu, delx, axis=3)
    dhvdy = np.gradient(hv, dely, axis=2)
    
    #dhudx = (hu[:,:,:-1,1:] - hu[:,:,:-1,:-1])/delx
    #dhvdy = (hv[:,:,1:,:-1] - hv[:,:,:-1,:-1])/dely
    divhu = dhudx + dhvdy
    divhu_tave = np.mean(divhu, axis=0)
    #h_tave = np.mean(c.cpd*T + L_v*QV, axis=0)
    #s_tave = np.mean(c.cpd*T, axis=0)
    #s_tave = s_tave + g*z3D
    #h_tave = h_tave + g*z3D
    h_tave = np.mean(h, axis=0)
    s_tave = np.mean(s, axis=0)
    
    #EDIT: this critical w threshold determines where convective region is. 
    #What is a good way to choose this value? (want to neglect gravity wave perturbations)
    Wcrit=0.01
    db=16
    totpoints=nx*ny
    
    rhoz = p/(c.Rd*np.mean(np.mean(T, axis=3), axis=2))
    
    rhoz_tave = np.mean(rhoz, axis=0)
    rhoz3D = np.zeros((ny, nx, nz))
    rhoz3D[:,:,:] = rhoz_tave
    rhoz3D = rhoz3D.T
    
    divhu_tave = np.multiply(divhu_tave, rhoz3D)
    
    sigma = np.zeros(z.size)
    l_ms = np.zeros(z.size)
    v_ns = np.zeros(z.size-1)
    divs = np.zeros(z.size-1)
    divh = np.zeros(z.size-1)
    div = np.zeros(z.size-1)
    Aup = np.zeros(z.size)
    h_convs = np.zeros(z.size-1)
    s_convs = np.zeros(z.size-1)
    
    
    

    
    #Wh_tave = np.multiply(Wtave, h_tave)
    #Ws_tave = np.multiply(Wtave, s_tave)
    
    #Wh_tave = blockave3D(Wh_tave, db)
    #Ws_tave = blockave3D(Ws_tave, db)
    
    #block averaging
    Wtave = blockave3D(Wtave, db)
    h_tave = blockave3D(h_tave, db)
    s_tave = blockave3D(s_tave, db)
    divhu_tave = blockave3D(divhu_tave, db)
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
        #Wh_down = Wh_tave[i,Wz > Wcrit]
        #Wh_up = Wh_tave[i+1, Wz > Wcrit]
        #Ws_down = Ws_tave[i, Wz > Wcrit]
        #Ws_up = Ws_tave[i+1, Wz > Wcrit]
        h_conv = np.mean(h_tave[i,Wz > Wcrit])
        s_conv = np.mean(s_tave[i, Wz > Wcrit])
        #vertical mass flux in convective region
        massflux = (rhoz_tave[i+1]*Wzup_crit - rhoz_tave[i]*Wz_crit)
        hflux = divhu_tave[i, Wz > Wcrit]
        #verthflux = (rhoz_tave[i+1]*Wh_up - rhoz_tave[i]*Wh_down)
        #vertsflux = (rhoz_tave[i+1]*Ws_up - rhoz_tave[i]*Ws_down)
        #divh[i] = rhoz_tave[i]*np.nansum(hflux)*((delx*dely)/Aup[i])
        divh[i] = np.nansum(hflux)*((delx*dely)/Aup[i])
        #divh[i] = rhoz_tave[i]*np.nanmean(hflux)
        #total vertical mass flux convergence in convective region
        #divs[i] = -np.sum(vertmassflux)*delx*dely*(db**2)
        div[i] = -np.nanmean(massflux)
        #divs[i] = -np.nanmean(vertsflux)
        #divh[i] = -np.nanmean(verthflux)
        
        #radial wind component (normal to convective region, assuming region is ciruclar)
        v_ns[i] = divs[i]/(2*np.pi*l_ms[i]*rhoz_tave[i])
        h_convs[i] = h_conv
        s_convs[i] = s_conv
        
    #sigmabar = np.mean(sigma)
    divbar = np.nanmean(np.abs(div))
    #l_mbar = np.mean(l_ms/1e3)
    #v_nbar = np.mean(v_ns[np.isfinite(v_ns)])
    divhbar = np.nanmean(np.abs(divh))
    
    BLi = np.where(p < 950*1e2)[0][0]
    if domsize == 768:
        p_outi = np.where(z > 13500)[0][0]
    elif domsize == 1536:
        p_outi = np.where(z > 13800)[0][0]
    elif domsize == 3072:
        p_outi = np.where(z > 14000)[0][0]
    else:
        p_outi = np.where(z > 13500)[0][0]
    MSE_BL = np.nanmean(h_convs[:BLi])
    MSE_t = np.mean(h_convs[p_outi-1:p_outi+1])
    s_BL = np.nanmean(s_convs[:BLi])
    s_t = np.mean(s_convs[p_outi-1:p_outi+1])

    delhtrop = MSE_BL - MSE_t
    delstrop = s_BL - s_t
    
    delz = np.diff(z)

        
    MSE_BL = np.nanmean(h_convs[:BLi])
    MSE_t = h_convs[p_outi]
    s_BL = np.nanmean(s_convs[:BLi])
    s_t = s_convs[p_outi]
    
    delhtrop = MSE_BL - MSE_t
    delstrop = s_BL - s_t
    
    #NGMS = np.nansum(divh*delz)/np.nansum(divs*delz)
        
    #print 'vertically averaged sigma', sigmabar
    #print 'vertically averaged l_m (km)', l_mbar
    #print 'vertically averaged v_n (m/s)', v_nbar
    #print 'vertically averaged mass flux (kg/s)', np.mean(divs)
    print 'vertically averaged divergence of MSE (W/m^2)', np.nanmean(divh)
    print 'vertical integral of MSE divergence (W/m^2)', np.nansum(divh*delz)
    #print 'vertical integral of DSE divergence (W/m^2)', np.nansum(divs*delz)
    #print 'Normalized GMS', NGMS
    print 'vertically averaged mass flux magnitude (kg/m^2)', divbar
    print 'delh_trop', delhtrop
    print 'MSE_BL', MSE_BL
    print 'MSE_t', MSE_t
    print 'dels_trop', delstrop
    print 's_BL', s_BL
    print 's_t', s_t
        
        
    #f, axarr = plt.subplots(3,1)
#    plt.figure(1)
#    ax = plt.gcf().gca()
#    if db == 1:
#        plt.title('convective fractional area (W$_{{c}}$ > {:2.3f} m/s)'.format(Wcrit))
#    else:
#        plt.title('convective fractional area (W$_{{c}}$ > {:2.3f} m/s), block-averaging over ({:2.0f} km)$^2$'.format(Wcrit, db*(np.diff(x)[0])/1e3))
#    plt.plot(sigma, z/1e3, '-x', color=colors[k], label='{:d} km, averaged over day {:2.1f} to {:2.1f}'.format(domsize, times[0], times[-1]))
#    ax.set_xlabel(r'$\sigma$')
#    plt.legend(loc='best')
#    plt.ylabel('z (km)')
#    #axarr[0,].plot(sigma, z/1e3)
#    #axarr[0,].set_xlabel('sigma')
#    #axarr[0,].set_ylabel('z (km)')
#    
#    ax.set_ylim(0, 20)
#    plt.savefig(fout + 'sigma_blkavg_day250_{:d}day_db{:d}Wc{:2.3f}.pdf'.format(nave, db, Wcrit))
#    
#    plt.figure(2)
#    ax = plt.gcf().gca()
#    if db == 1:
#        plt.title('normal velocity along convective region edge, W$_{{c}}$ = {:4.2} m/s'.format(Wcrit))
#    else:
#        plt.title('normal velocity along convective region edge, W$_{{c}}$ = {:2.1f} m/s, block-averaging over ({:2.0f} km)$^2$'.format(Wcrit, db*(np.diff(x)[0])/1e3))
#    plt.plot(v_ns, z[:-1]/1e3, '-x', color=colors[k], label='{:d} km, averaged over day {:2.1f} to {:2.1f}'.format(domsize, times[0], times[-1]))
#    plt.xlabel('velocity (m/s)')
#    plt.ylabel('z (km)')
#    plt.legend(loc='best')
#    ax.set_ylim(0, 20)
#    #ax.set_xlim(-10, 300)
#    
#    plt.savefig(fout + 'vnormal_blkavg_day250_{:d}day_db{:d}Wc{:2.3f}.pdf'.format(nave, db, Wcrit))
#
#    plt.figure(3)
#    ax = plt.gcf().gca()
#    if db == 1:
#        plt.title(r'horizontal mass flux out of convective region, W$_{{c}}$ = {:4.2f} m/s'.format(Wcrit))
#    else:
#        plt.title(r'horizontal mass flux out of convective region, W$_{{c}}$ = {:4.2f} m/s, block-averaging over ({:2.0f} km)$^2$'.format(Wcrit, db*(np.diff(x)[0])/1e3))
#    plt.plot(divs/divbar, z[:-1]/1e3, '-x', color=colors[k],  label='{:d} km, averaged over day {:2.1f} to {:2.1f}'.format(domsize, times[0], times[-1]))
#    plt.xlabel('mass flux relative to vertically-averaged magnitude of mass flux')
#    plt.ylabel('z (km)')
#    plt.legend(loc='best')
#    ax.set_ylim(0, 20)
    
    plt.savefig(fout + 'massflux_blkavg_day250_{:d}day_db{:d}Wc{:2.3f}.pdf'.format(nave, db, Wcrit))
    
    plt.figure(6)
    ax = plt.gcf().gca()
    if db == 1:
        plt.title(r'horizontal MSE flux out of convective region, $w_{{c}}$ = {:4.2f} m/s'.format(Wcrit))
    else:
        plt.title(r'horizontal MSE flux out of convective region, $w_{{c}}$ = {:4.2f} m/s, block-averaging over ({:2.0f} km)$^2$'.format(Wcrit, db*(np.diff(x)[0])/1e3))
    plt.plot(divh/divhbar, z[:-1]/1e3, '-x', color=colors[k],  label='{:d} km, averaged over day {:2.0f} to {:2.0f}'.format(domsize, times[0], times[-1]))
    plt.xlabel(r'$\frac{(\nabla \cdot \rho{\bf u}h)_{conv}}{\overline{|\nabla \cdot \rho{\bf u}h|}_{conv}}$')
    plt.ylabel('z (km)')
    plt.legend(loc='best')
    ax.set_ylim(0, 20)
    
    plt.savefig(fout + 'MSEflux_blkavg_day250_{:d}day_db{:d}Wc{:2.3f}.pdf'.format(nave, db, Wcrit))
    
    plt.figure(7)
    ax = plt.gcf().gca()
    if db == 1:
        plt.title(r'MSE in convective region, $w_{{c}}$ = {:4.2f} m/s'.format(Wcrit))
    else:
        plt.title(r'MSE in convective region, $w_{{c}}$ = {:4.2f} m/s, block-averaging over ({:2.0f} km)$^2$'.format(Wcrit, db*(np.diff(x)[0])/1e3))
    plt.plot(h_convs, z[:-1]/1e3, '-x', color=colors[k],  label='{:d} km, averaged over day {:2.0f} to {:2.0f}'.format(domsize, times[0], times[-1]))
    plt.xlabel('MSE (J)')
    plt.ylabel('z (km)')
    plt.legend(loc='best')
    ax.set_ylim(0, 17)
    
    plt.savefig(fout + 'MSE_blkavg_day250_{:d}day_db{:d}Wc{:2.3f}.pdf'.format(nave, db, Wcrit))
    
    plt.figure(8)
    ax = plt.gcf().gca()
    if db == 1:
        plt.title(r'DSE in convective region, $w_{{c}}$ = {:4.2f} m/s'.format(Wcrit))
    else:
        plt.title(r'DSE in convective region, $w_{{c}}$ = {:4.2f} m/s, block-averaging over ({:2.0f} km)$^2$'.format(Wcrit, db*(np.diff(x)[0])/1e3))
    plt.plot(s_convs, z[:-1]/1e3, '-x', color=colors[k],  label='{:d} km, averaged over day {:2.0f} to {:2.0f}'.format(domsize, times[0], times[-1]))
    plt.xlabel('MSE (J)')
    plt.ylabel('z (km)')
    plt.legend(loc='best')
    ax.set_ylim(0, 17)
    
    plt.savefig(fout + 'MSE_blkavg_day250_{:d}day_db{:d}Wc{:2.3f}.pdf'.format(nave, db, Wcrit))
    
    #vertically smooth the divergence
    
    #def moving_average(a, n=3):
    #    ret = np.cumsum(a, dtype=float)
    #    ret[n:] = ret[n:] - ret[:-n]
    #    return ret[n - 1:] / n
    #
    #n=5
    #divsmooth = moving_average(divs/divbar, n)
    #
    #plt.figure(4)
    #ax = plt.gcf().gca()
    #
    #trunc = (n-1)/2
    #if db == 1:
    #    plt.title(r'vertically smoothed horizontal mass flux out of convective region, W$_{{c}}$ = {:4.2f} m/s, window size = {:d}'.format(Wcrit, n))
    #else:
    #    plt.title(r'vertically smoothed horizontal mass flux out of convective region, W$_{{c}}$ = {:4.2f} m/s, window size = {:d}, block-averaging over ({:2.0f} km)$^2$'.format(Wcrit, n, db*(np.diff(x)[0])/1e3))
    #plt.plot(divsmooth, z[trunc:-(trunc+1)]/1e3, '-x', color=colors[k], label='{:d} km, averaged over day {:2.1f} to {:2.1f}'.format(domsize, times[0], times[-1]))
    #plt.xlabel('mass flux relative to vertically-averaged magnitude of mass flux')
    #plt.ylabel('z (km)')
    #plt.legend(loc='best')
    #ax.set_ylim(0, 20)
    #
    #plt.savefig(fout + 'smoothmassflux_blkavg_day250_{:d}day_db{:d}Wc{:2.3f}.pdf'.format(nave, db, Wcrit))
    
#plt.figure(1)
#plt.close()
#
#plt.figure(2)
#plt.close()

plt.figure(6)
plt.close()

plt.figure(7)
plt.close()

plt.figure(8)
plt.close()
        
        
    


#axarr[1,].plot(l_ms/1e3, z/1e3)
#axarr[1,].set_xlabel('radius of convective region (km)')
#axarr[1,].set_ylabel('z (km)')

#axarr[2,].plot(Aup/(1e6), z/1e3)
#axarr[2,].set_xlabel(r'area of convective region (km$^2$)')
#axarr[2,].set_ylabel('z (km)')







    
    
    
    
    

    
    
    



