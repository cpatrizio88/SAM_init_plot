from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
import SAM_init_plot.block_fns
import SAM_init_plot.misc_fns
import thermolib
from SAM_init_plot.misc_fns import radprof3D
from SAM_init_plot.block_fns import blockave2D, blockave3D, blockxysort2D
from thermolib.constants import constants
from thermolib.wsat import wsat
import site
import sys
#site.addsitedir('\\Users\\Casey\\Dropbox\\research\\code\\thermlib')
site.addsitedir('/Users/cpatrizio/repos/thermolib/')
import thermo
from findTmoist import findTmoist


c=constants()

matplotlib.rcParams.update({'font.size': 26})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 24})
matplotlib.rcParams.update({'mathtext.fontset': 'cm'})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fpathSTAT = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
fout='/Users/cpatrizio/Google Drive/MS/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_RADIAL/'

nc_in3D1 = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]
nc_in2D1 = glob.glob(fpath2D + '*256x256*3000m*day230to250*302K.nc')[0]
nc_inSTAT1 = glob.glob(fpathSTAT + '*256x256*3000m*250days*302K.nc')[0]

nc_in2D2 = glob.glob(fpath2D + '*512x512*3000m*day180to195*302K.nc')[0]
nc_in3D2 = glob.glob(fpath3D + '*512x512*3000m*day180to195*302K.nc')[0]
nc_inSTAT2 = glob.glob(fpathSTAT + '*512x512*3000m*195days*302K.nc')[0]

nc_in3D3 = glob.glob(fpath3D + '*1024x1024*day170to180*302K.nc')[0]
nc_in2D3 = glob.glob(fpath2D + '*1024x1024*3000m*day170to180*302K.nc')[0]
nc_inSTAT3 = glob.glob(fpathSTAT + '*1024x1024*3000m*220days*302K.nc')[0]

domsizes=[768, 1536, 3072]
colors=['k', 'r', 'g']


nc_STATs = [nc_inSTAT1, nc_inSTAT2, nc_inSTAT3]
nc_2Ds = [nc_in2D1, nc_in2D2, nc_in2D3]
nc_3Ds = [nc_in3D1, nc_in3D2, nc_in3D3]

for i, domsize in enumerate(domsizes):
    if domsize == 768:
        nave=8
    elif domsize == 1536:
        nave=8
    else: 
        nave=8

    ntave2D=24
    ntave3D=4
    t2 = -1
    t3 = -1

    aveperiod2D = nave*ntave2D
    aveperiod3D = nave*ntave3D

    print 'loading 2D variables'
    nc_data2D = Dataset(nc_2Ds[i])
    varis2D = nc_data2D.variables
    
    PW = varis2D['PW'][t2-aveperiod2D:t2,:,:]
    
    print 'loading 3D variables'
    nc_data3D = Dataset(nc_3Ds[i])
    varis3D = nc_data3D.variables
    
    print 'loading STAT variables'
    nc_dataSTAT = Dataset(nc_STATs[i])
    varisSTAT = nc_dataSTAT.variables

    theta = varisSTAT['THETA'][t2-aveperiod2D:t2,:]
    RELH = varisSTAT['RELH'][t2-aveperiod2D:t2,:]
    RELH_tave = np.mean(RELH, axis=0)
    theta_tave = np.mean(theta, axis=0)

    t3D = varis3D['time'][t3-aveperiod3D:t3]
    t2D = varis2D['time'][t2-aveperiod2D:t2]
    x = varis3D['x'][:]
    y = varis3D['y'][:]
    z = varis3D['z'][:]
    p = varis3D['p'][:]
    p = p*1e2

    #varnames = ['QRAD', 'W', 'U', 'QN', 'QV']
    #varnames=['QV', 'TABS']
    varnames=['QV']
    varnames = ['TABS']

    #calculate relative humidity 
    
    #averaging time period for 2D & 3D fields
    ntave2D=24
    ntave3D=4
    
    nt3D = t3D.size
    nt2D = t2D.size
    nz = z.size
    nx = x.size
    ny = y.size
    
    xx, yy = np.meshgrid(x, y)
    times = np.arange(t3D[0], np.max(t3D))

    #2D fields
    #ntrunc = PW.shape[0]%ntave2D
    #PW = PW[ntrunc:,:,:]
    #PW_tmp = PW.reshape(nt2D/ntave2D, ntave2D, nx, ny)
    #PW_tave = np.mean(PW_tmp, axis=1)
    
    #time period to look at
    #t=-25*ntave3D
    #nave=5
    #nave2D=5*ntave2D
    #nave3D=5*ntave3D
    #number of blocks to find PW max
    db=16
    print 'calulating max PW, blocked'
    #PW_t = np.mean(PW_tave[t-nave:t,:,:], axis=0)
    #PW_t = np.mean(PW[t-nave2D:t,:,:], axis=0)
    PW_tave = np.mean(PW, axis=0)
    
    PW_blocked = blockave2D(PW_tave, db)
    PWxy_sorted = blockxysort2D(PW_tave, xx, yy, db)
    
    PWsort = PWxy_sorted.keys()
    PWxycoords = PWxy_sorted.values()
    
    mcenter = PWxycoords[-1]
    binwidth=5e3
    rbins = np.arange(0, domsize*1e3, binwidth)
    
    #FIND BL TOP
    
    dthetadz = (theta_tave[1:]-theta_tave[:-1])/np.diff(z)
    dRHdz = (RELH_tave[1:] - RELH_tave[:-1])/np.diff(z)
    
    dthetadz_crit = 1e-3
    
    BLi = np.where(dthetadz > dthetadz_crit)[0][0]
    BLi = np.where(z > 500)[0][0]
    
    p_BL = p[BLi]
    z_BL = z[BLi]
    zBL = z[:BLi]
    
    nbins=[rbins, zBL]
    
    print 'calculating z-r plots'
    for varname in varnames:
        print 'loading', varname
        vari = varis3D[varname]
        field = varis3D[varname][t3-aveperiod3D:t3,:,:,:]
        field_tave = np.mean(field, axis=0)
        #calculate relative humidity 
        if varname == 'TABS':
            varname = 'RH'
            QV = varis3D['QV'][t3-aveperiod3D:t3,:,:,:]
            QV_tave = np.mean(QV, axis=0)
            RH_tave = np.zeros(QV_tave.shape)
            for j, plev in enumerate(p):
                wsat_tave = 1000*wsat(field_tave[j,:,:], plev) #convert to g/kg
                RH_tave[j,:,:] = 100*(QV_tave[j,:,:]/wsat_tave) #convert to percent
            field_tave = RH_tave
        #calculate wind speed
        if varname == 'U':
            V = varis3D['V'][t3-aveperiod3D:t3,:,:,:]
            #V = V[ntrunc:,:,:,:]
            #U = varis3D['U'][t3-aveperiod3D:t3,:,:,:]
            #U = U[ntrunc:,:,:,:]
            field = np.sqrt(np.power(V, 2) + np.power(field, 2))
            #field_tmp = field.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
            #field_tave = np.mean(field_tmp, axis=1)
            print 'calculating temporal averages'
            #average over nave days 
            #field_t = np.mean(field[t-nave3D:t,:,:,:], axis=0)
            #calculate adiabatic warming + QRAD balance 
        if varname == 'QRAD':
            Qr = field_tave
            varname = 'ADBWARMplusQRAD'
            w = varis3D['W'][t3-aveperiod3D:t3,:,:,:]
            T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
            #w_t = np.mean(w[t-nave3D:t,:,:,:], axis=0)
            w_tave = np.mean(w, axis=0)
            #T_t = np.mean(T[t-nave3D:t,:,:,:], axis=0)
            T_tave = np.mean(T, axis=0)
            #z3D = np.zeros(field_t.shape)
            z3D = np.zeros(field_tave.shape)
            z3D = z3D.T
            z3D[:,:,:] = z
            z3D = z3D.T
            delz3D = np.zeros((nz-1, nx, ny))
            delz3D = delz3D.T
            delz3D[:,:,:] = np.diff(z)
            delz3D = delz3D.T
            s = c.cp*T_tave + c.g*z3D
            dsdz = (s[1:,:,:] - s[:-1,:,:])/delz3D
            adbwarm = -(1./c.cp)*np.multiply(w_tave[:-1], dsdz)
            adbwarm = adbwarm*(3600*24) #convert from K/s to K/day
            field_tave = adbwarm + Qr[:-1,:,:]
            znew = z[:-1]
        
        
        #fieldmeans has shape (redges.size, zedges.size) = nbins
        print '2d contouring'
        if varname == 'ADBWARMplusQRAD':
            redges, zedges, fieldmeans = radprof3D(field_tave, xx, yy, znew, mcenter, nbins=nbins)
        else:
            redges, zedges, fieldmeans = radprof3D(field_tave, xx, yy, z, mcenter, nbins=nbins)
        rbin_centers = (redges[1:] + redges[:-1])/2.
        zbin_centers = (zedges[1:] + zedges[:-1])/2.
        
        #rr has shape (zedges.size, redges.size) = nbins.T
        rr, zz = np.meshgrid(rbin_centers, zbin_centers)
        
        fieldmeans = np.transpose(fieldmeans)
        
        if varname == 'RH':
            units='%'
        #    cb.set_label(units)
        else:
            units = vari.units.strip()
        #    cb.set_label(units)
        
        
        
        #plt.figure()
        #ax=plt.gcf().gca()
        #vmin=np.ma.masked_invalid(fieldmeans).min()
        #vmax=np.ma.masked_invalid(fieldmeans).max()
        ##EDIT: COLORBAR FORMATTING. IF NOT MONOTONIC, SWITCH TO DIVERGING COLORBAR, OTHERWISE
        ##USE SEQUENTIAL COLORBAR. SET COLOR MAPPING LIMITS MANUALLY. 
        #if (vmin < 0 and vmax > 0):
        #   if np.abs(vmin) < np.abs(vmax):
        #      vmin = vmin - (vmax - np.abs(vmin))
        #   else:
        #      vmax = vmax + (np.abs(vmin) - vmax)
        #   cmap = cm.RdBu_r
        #else:
        #   cmap = cm.YlGnBu
        #if varname == 'W':
        #   vmax=.4
        #   vmin=-.4
        #   cmap = cm.RdBu_r
        #if varname == 'ADBWARMplusQRAD':
        #   vmin=-6
        #   vmax=6
        #   cmap=cm.RdBu_r
        #if varname == 'U':
        #   vmin=0
        #   vmax=9
        #if varname == 'QN':
        #   vmin=0
        #   vmax=0.2
        #print 'plotting'
        #fieldmeans = np.transpose(fieldmeans)
        ##extent=[redges[0]/(domsize*1e3), redges[-1]/(domsize*1e3), zedges[0]/1e3, zedges[-1]/1e3]
        ##plt.contourf(rr/(1e3*domsize), zz/(1e3), np.transpose(fieldmeans), cmap=cm.RdYlBu_r)
        ##plt.imshow(np.transpose(fieldmeans), origin='lower', interpolation='nearest', aspect='auto', extent=extent, cmap=cm.RdYlBu_r) 
        #if varname == 'QV':
        #     plt.pcolormesh(rr/(1e3*domsize), zz/(1e3), fieldmeans, vmin=vmin, vmax=vmax, cmap=cmap, norm=matplotlib.colors.LogNorm())
        #else:
        #      plt.pcolormesh(rr/(1e3*domsize), zz/(1e3), fieldmeans, vmin=vmin, vmax=vmax, cmap=cmap)
        #plt.xlabel('fractional distance from moist region center, relative to domain size')
        #plt.ylabel('z (km)')
        #cb=plt.colorbar()
        #cb.set_label(units)
        ##ax.set_ylim(0,26)
        ##ax.set_ylim(0, z_BL/1e3)
        #ax.set_xlim(0, 0.7)
        #plt.title('{:s} ({:s}), day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$, x bin width = {:2.2f} km'.format(varname, units, t3D[0], t3D[-1], domsize, binwidth/1e3))
        ##plt.title('{:s} ({:s}), day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$, x bin width = {:2.2f}'.format(varname, vari.units.strip(), t3D[t-nave*ntave3D], t3D[t], domsize, 1./nbins[0]))
        #plt.savefig(fout + '{:s}radialxsection_day{:3.0f}to{:3.0f}.pdf'.format(varname, t3D[0], t3D[-1]))
        ##plt.savefig(fout + '{:s}radialxsection_day{:3.0f}to{:3.0f}.pdf'.format(varname, t3D[t-nave*ntave3D], t3D[t]))
        #plt.close()
        g = 9.81
        sig = 5.67e-8
        T_s = 302 
        L_v = 2.257*1e6
        
        S0 = 413.98 #solar constant
        theta = 50.5 #zenith angle 
        S = np.cos((theta*np.pi/180))*S0 
        
        c_E = 0.001
        rho_w = 1000 #density of water
        
        
        d = 100000e3 #
        
        p_s = 1000e2 #surface pressure (Pa)
        p_t = 200e2 #tropopause (Pa)
        p_BL = 900e2 #boundary layer top (Pa)
        T_s = 302
        delz_BL = z[BLi]
        c_E = 0.001
        l_d = domsize*1e3*(1./np.sqrt(2))
        #r = np.linspace(0, l_d, 1e6)
        q_sat = wsat(T_s, p_s) #mixing ratio above sea surface (100% saturated)    
    
        q_sat = wsat(T_s, p_s) #mixing ratio above sea surface (100% saturated)    
        thetae0 = thermo.theta_e(T_s, p_s, q_sat, 0) #theta_e in moist region 
        T_BLtop = findTmoist(thetae0, p_BL) #temperature of boundary layer top 
        T_t = findTmoist(thetae0, p_t) #temperature of tropopause (outflow region)
        T_BL = (T_s + T_BLtop)/2. #temperature of boundary layer, consistent with well-mixed assumption (linear mixing)
        q_FA = wsat(T_t, p_t) #free troposphere water vapor mixing ratio
        
        zhat = delz_BL/c_E
    
    
        q_BL = q_sat + 2*zhat*(q_FA - q_sat)/(l_d**2 - rbin_centers**2)*(rbin_centers + zhat - (zhat + l_d)*np.exp((rbin_centers - l_d)/zhat))
        
        RH_c = q_BL/q_sat
        
        if varname == 'QV':
           titlename = r'$q_v$'
        if varname == 'RH':
           titlename = 'RH'

    
    
        #plot vertical means
        plt.figure(1)
        fieldmeanbar = np.mean(fieldmeans, axis=0)
        fieldaverage = np.nanmean(fieldmeanbar)
        
        plt.plot(rbin_centers/(1e3*domsize), fieldmeanbar, color=colors[i], label='{:d} km, day {:3.0f} to {:3.0f}'.format(domsize, t3D[0], t3D[-1]))
        plt.plot(rbin_centers/(1e3*domsize), 100*RH_c, color=colors[i], alpha=0.5, linewidth=2, label = '{:d} km, simple model'.format(domsize))
        plt.xlabel(r'$\hat{r}$')
        plt.ylim(0, 100)
        plt.xlim(0, 1/np.sqrt(2))
        #plt.xlabel('fractional distance from moist region center, relative to domain size')
        plt.ylabel('{:s} ({:s})'.format(titlename, units))
        plt.title('BL {:s}' .format(titlename))
        plt.savefig(fout + '{:s}radialvertmean_day{:3.0f}to{:3.0f}.pdf'.format(varname, t3D[0], t3D[-1]))
        plt.legend()
        plt.show()
        
        #qBLbar_s =  17.5 
        
        if varname == 'QV':
            print 'simple model average qBL = {:3.1f} g/kg'.format(np.mean(q_BL)*1e3)
            print 'simple model average RH = {:3.3f} %'.format(RH_c)
            print 'simulation average {:s} = {:3.1f} g/kg'.format(varname, fieldaverage)
            

        
    
    
    
    
    


