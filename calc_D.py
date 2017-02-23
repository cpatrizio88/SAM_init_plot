from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from SAM_init_plot.block_fns import blockave2D, blockave3D
from thermolib.constants import constants
import matplotlib.colors as colors

c = constants()

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fpathSTAT = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'


nc_inSTAT1 = glob.glob(fpathSTAT + '*256x256*3000m*250days*302K.nc')[0]
nc_in2D1 = glob.glob(fpath2D + '*256x256*3000m*day230to250*302K.nc')[0]
nc_in3D1 = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]


nc_inSTAT2 = glob.glob(fpathSTAT + '*512x512*3000m*180days*302K.nc')[0]
nc_in2D2 = glob.glob(fpath2D + '*512x512*3000m*day180to195*302K.nc')[0]
nc_in3D2 = glob.glob(fpath3D + '*512x512*3000m*day180to195*302K.nc')[0]

nc_inSTAT3 = glob.glob(fpathSTAT + '*1024x1024*3000m*180days*302K.nc')[0]
nc_in2D3 = glob.glob(fpath2D + '*1024x1024*3000m*day170to180*302K.nc')[0]
nc_in3D3 = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]

fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STAT/'


domsizes = [768, 1536, 3072]
nc_STATs = [nc_inSTAT1, nc_inSTAT2, nc_inSTAT3]
nc_2Ds = [nc_in2D1, nc_in2D2, nc_in2D3]
nc_3Ds = [nc_in3D1, nc_in3D2, nc_in3D3]

#nc_3Ds = [nc_in3D1]

pcolors = ['k', 'r', 'g']


tstarts = [0,90,90]

for i, nc_f3D in enumerate(nc_3Ds):
        
    print 'domsize', domsizes[i]
    
    print 'loading netCDF files'
    
    db=1
    W_c = 0.5
    
    times = np.array([])
    
    #varname='W'
    
    convfracs = np.array([])

    #for nc_in3D in nc_f3D:
    #nc_data3D = Dataset(nc_in3D) 
        
    print 'loading', nc_f3D
    
    nc_data3D = Dataset(nc_f3D)
    varis3D = nc_data3D.variables
    
    nc_data2D = Dataset(nc_2Ds[i])
    varis2D = nc_data2D.variables
    
    nc_dataSTAT = Dataset(nc_STATs[i])
    varisSTAT = nc_dataSTAT.variables
    
    P = varis2D['Prec'][:]
    
    rho_w = 1000
    #convert from mm/day to kg m^-2 s^-1
    P = P*(rho_w*.001)/(86400)
    
    Pbar = np.mean(np.mean(P, axis=1), axis=1)
    
    units = r'kg m$^{-2}$ s$^{-1}$'
    
    x = varis3D['x'][:]
    y = varis3D['y'][:]
    z = varisSTAT['z'][:]
    p = varisSTAT['p'][:]
    p=p*1e2
    z_BL = 2e3
    
    BLi = np.where(z > z_BL)[0][0]
    
    
    #nt3D = t3D.size
    nz = z.size
    nx = x.size
    ny = y.size
    
    nxprime = nx / db
    nyprime = ny / db
    
    #averaging time period for 2D & 3D fields
    ntave2D=24
    ntave3D=4

    t3D = varis3D['time'][:]
    nt = len(t3D)
    
    t2D = varis2D['time'][:]
    #wfield = varis3D[varname][:]
    
    P2D = varis2D['Prec'][:]
    
    #ntrunc = nt2D%ntave2D
    #ntrunc=0
    #field = field[ntrunc:,:,:]
    #units = varis3D[varname].units
    convfrac = np.zeros(nt)
    qvflux = np.zeros(nt)
    mflux = np.zeros(nt)
    qvBL = np.zeros(nt)
    d = np.zeros(nt)
    Pconv = np.zeros(nt)

    for k, t in enumerate(t3D):
        k2 = (ntave2D/ntave3D)*k
        #Pfield_t = P2D[k2,:,:]
        #wfield_t = wfield[k,:,:,:]
        w = varis3D['W'][k,:,:,:]
        T = varis3D['TABS'][k,:,:,:]
        QV = varis3D['QV'][k,:,:,:]
        QV = QV*1e-3
        nt = T.shape[0]
        p3D = np.zeros((ny, nx, nz))
        p3D[:,:,:] = p
        p3D = p3D.T
        rho = p3D/(T*c.Rd)
        massflux = np.multiply(w, rho)
        massflux_blockave = blockave3D(massflux, db)
        wfield_blockave = blockave3D(w, db)
        QV_blockave = blockave3D(QV, db)
        P_blockave = blockave2D(P[k2,:,:], db)
        #P_blockave = blockave2D(Pfield_t, db)
        massflux_BLave = np.mean(massflux_blockave[:BLi,:,:], axis=0) 
        QV_BLave = np.mean(QV_blockave[:BLi,:,:], axis=0)
        w_BLave = np.mean(wfield_blockave[:BLi,:,:], axis=0)
        
        convBL_points = w_BLave > W_c
        
        massflux_convave = np.mean(massflux_BLave[convBL_points])
        QV_convave = np.mean(QV_BLave[convBL_points])
        P_convave = np.mean(P_blockave[convBL_points])

        qvflux[k] = massflux_convave*QV_convave
        mflux[k] = massflux_convave
        qvBL[k] = QV_convave
        Pconv[k] = P_convave
        d[k] = qvflux[k] - Pconv[k]
        
        #plt.close()
        convfrac[k] = len(w_BLave[convBL_points])/(1.*nxprime*nyprime)
        #subsfrac[i] = 1 - convfrac[i]
    
    #        buoyflux = np.vstack((buoyflux, buoyflux_temp)) if buoyflux.size else buoyflux_temp
    convfracs = np.concatenate((convfracs, convfrac)) if convfracs.size else convfrac
    times = np.concatenate((times, t3D)) if times.size else t3D
    
    #ntimes = convfracs.shape[0]
    
    print 'plotting'
    
    #calculate 2D (CRH rank, time) map of CRH frequency  ###
    #times = np.arange(tstarts[i], tstarts[i]+ntimes)

    w = 4*5
    w2 = 24*5
    
    convfracs_smooth = moving_average(convfracs, w)
    qvflux_smooth = moving_average(qvflux, w)
    mflux_smooth = moving_average(mflux, w)
    d_smooth = moving_average(d, w)
    Pbar_smooth = moving_average(Pbar, w2)
    Pconv_smooth = moving_average(Pconv, w)
    
    plt.figure(1)
    plt.plot(times, convfracs, color=pcolors[i], alpha=0.7, label='{:d} km'.format(domsizes[i]))
    plt.plot(times[w/2:-w/2+1], convfracs_smooth, color=pcolors[i], linewidth=1.5)
    #plt.plot(times, 1-convfracs, 'k--')
    plt.xlabel('time (days)')
    plt.ylabel('convective fractional area')
    #plt.ylim((0, 0.10))
    if db == 1:
        plt.title('BL convective fractional area (W  > W$_c$ = {:3.2f} m/s)'.format(W_c))
    else:
        plt.title('BL convective fractional area (W > W$_c$ = {:3.2f} m/s), average over ({:2.0f} km)$^2$ blocks '.format(W_c, db*np.diff(x)[0]/1e3))
    plt.savefig(fout + 'fracconvNEW250days_db{:d}.pdf'.format(db))
    
    plt.figure(2)
    plt.plot(times, qvflux, color=pcolors[i], alpha=0.7, label='{:d} km'.format(domsizes[i]))
    plt.plot(times[w/2:-w/2+1], qvflux_smooth, color=pcolors[i], linewidth=1.5)
    #plt.plot(times[w/2:-w/2+1], convfracs_smooth, color=pcolors[i], linewidth=1.5)
    #plt.plot(times, 1-convfracs, 'k--')
    plt.xlabel('time (days)')
    plt.ylabel(r'q$_v$ flux (kg m$^{-2}$ s$^{-1}$)')
    plt.gca().set_ylim(bottom=0)
    #plt.ylim((0, 0.10))
    if db == 1:
        plt.title(r'vertical q$_v$ flux in convective region BL (W  > W$_c$ = {:3.2f} m/s)'.format(W_c))
    else:
        plt.title(r'vertical q$_v$ flux in convective region BL (W > W$_c$ = {:3.2f} m/s), average over ({:2.0f} km)$^2$ blocks '.format(W_c, db*np.diff(x)[0]/1e3))
    plt.savefig(fout + 'qvflux250days_db{:d}.pdf'.format(db))
    
    plt.figure(3)
    plt.plot(times, mflux, color=pcolors[i], alpha=0.7, label='{:d} km'.format(domsizes[i]))
    plt.plot(times[w/2:-w/2+1], mflux_smooth, color=pcolors[i], linewidth=1.5)
    #plt.plot(times[w/2:-w/2+1], convfracs_smooth, color=pcolors[i], linewidth=1.5)
    #plt.plot(times, 1-convfracs, 'k--')
    plt.xlabel('time (days)')
    plt.ylabel(r' $\rho w$ (kg m$^{-2}$ s$^{-1}$)')
    #plt.ylim((0, 0.10))
    if db == 1:
        plt.title(r'$\rho w$ in convective region BL (W  > W$_c$ = {:3.2f} m/s)'.format(W_c))
    else:
        plt.title(r'$\rho w$ in convective region BL (W > W$_c$ = {:3.2f} m/s), average over ({:2.0f} km)$^2$ blocks '.format(W_c, db*np.diff(x)[0]/1e3))
    plt.savefig(fout + 'massflux250days_db{:d}.pdf'.format(db))
    
    plt.figure(4)
    plt.plot(times, d, color=pcolors[i], alpha=0.7, label='{:d} km'.format(domsizes[i]))
    plt.plot(times[w/2:-w/2+1], d_smooth, color=pcolors[i], linewidth=1.5)
    plt.axhline(0, alpha=0.3)
    #plt.plot(times[w/2:-w/2+1], convfracs_smooth, color=pcolors[i], linewidth=1.5)
    #plt.plot(times, 1-convfracs, 'k--')
    plt.xlabel('time (days)')
    plt.ylabel(r' $d$ (kg m$^{-2}$ s$^{-1}$)')
    #plt.ylim((0, 0.10))
    if db == 1:
        plt.title(r'convective region detrainment parameter, $d$')
    else:
        plt.title(r'convective region detrainment parameter, $d$, fields averaged over ({:2.0f} km)$^2$ blocks'.format(db*np.diff(x)[0]/1e3))
    plt.savefig(fout + 'd250days_db{:d}.pdf'.format(db))
    
    plt.figure(5)
    plt.plot(t2D, Pbar, '--', color=pcolors[i], alpha=0.7, label='{:d} km, domain mean'.format(domsizes[i]))
    plt.plot(times, Pconv, color=pcolors[i], alpha=0.7, label = '{:d} km, convective'.format(domsizes[i]))
    plt.plot(times[w/2:-w/2+1], Pconv_smooth, color=pcolors[i], linewidth=1.5)
    #plt.plot(times[w/2:-w/2+1], convfracs_smooth, color=pcolors[i], linewidth=1.5)
    #plt.plot(times, 1-convfracs, 'k--')
    plt.xlabel('time (days)')
    plt.ylabel(r' $P$ (kg m$^{-2}$ s$^{-1}$)')
    plt.gca().set_ylim(bottom=0)
    plt.title('mean precipitation, P')
    plt.savefig(fout + 'P250days_db{:d}.pdf'.format(db))
    

  

plt.figure(1)
plt.legend(loc='best')
plt.savefig(fout + 'fracconv250days_db{:d}.pdf'.format(db))
plt.close()

plt.figure(2)
plt.legend(loc='best')
plt.savefig(fout + 'qvflux250days_db{:d}.pdf'.format(db))
plt.close()

plt.figure(3)
plt.legend(loc='best')
plt.savefig(fout + 'massflux250days_db{:d}.pdf'.format(db))
plt.close()

plt.figure(4)
plt.legend(loc='best')
plt.savefig(fout + 'd250days_db{:d}.pdf'.format(db))
plt.close()

plt.figure(5)
plt.legend(loc='best')
plt.savefig(fout + 'P250days_db{:d}.pdf'.format(db))
plt.close()

        

    