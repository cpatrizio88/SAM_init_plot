from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from thermolib.wsat import wsat
from thermolib.constants import constants
from SAM_init_plot.block_fns import blockave3D
import SAM_init_plot.misc_fns 
import gc

c = constants()

matplotlib.rcParams.update({'font.size': 30})
matplotlib.rcParams.update({'figure.figsize': (24, 18)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 22})

plt.style.use('seaborn-white')

fpath2D =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpathSTAT = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'

foutdata = '/Users/cpatrizio/data/SST302/'

#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to140_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to130_64vert_ubarzero_MOISTDRYPROFS/'

fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_MOISTDRYPROFS/'

nc_inSTAT1 = glob.glob(fpathSTAT + '*256x256*3000m*250days*302K.nc')[0]
nc_in2D1 = glob.glob(fpath2D + '*256x256*3000m*day230to250*302K.nc')[0]
nc_in3D1 = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]


nc_inSTAT2 = glob.glob(fpathSTAT + '*512x512*3000m*180days*302K.nc')[0]
nc_in2D2 = glob.glob(fpath2D + '*512x512*3000m*day180to195*302K.nc')[0]
nc_in3D2 = glob.glob(fpath3D + '*512x512*3000m*day180to195*302K.nc')[0]

nc_inSTAT3 = glob.glob(fpathSTAT + '*1024x1024*3000m*180days*302K.nc')[0]
nc_in2D3 = glob.glob(fpath2D + '*1024x1024*3000m*day170to180*302K.nc')[0]
nc_in3D3 = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]


domsizes = [768, 1536, 3072]
nc_STATs = [nc_inSTAT1, nc_inSTAT2, nc_inSTAT3]
nc_2Ds = [nc_in2D1, nc_in2D2, nc_in2D3]
nc_3Ds = [nc_in3D1, nc_in3D2, nc_in3D3]

#domsizes = [3072]
#nc_STATs = [nc_inSTAT3]
#nc_2Ds = [nc_in2D3]
#nc_3Ds = [nc_in3D3]

#domsizes = [768]
#nc_STATs = [nc_inSTAT1]
#nc_2Ds = [nc_in2D1]
#nc_3Ds = [nc_in3D1]

#domsizes = [1536]
#nc_STATS = [nc_inSTAT2]
#nc_2Ds = [nc_in2D2]
#nc_3Ds = [nc_in3D2]

#domsizes = [768, 1536]
#nc_STATs = [nc_inSTAT1, nc_inSTAT2]
#nc_2Ds = [nc_in2D1, nc_in2D2]
#nc_3Ds = [nc_in3D1, nc_in3D2]

colors = ['k', 'r', 'g']
#colors = ['k', 'r']
#colors = ['k', 'r']

nave=10
ntave2D=24
ntave3D=4
#t2 = -35*ntave2D
#t3 = -35*ntave3D

t2s = [-1, -1, -1]
t3s = [-1, -1, -1]

aveperiod2D = nave*ntave2D
aveperiod3D = nave*ntave3D

varnames = ['W']

#varnames = ['CLOUD AMOUNT']

#varnames = ['QV']

#varnames = ['URADIAL']

#varnames = ['TABS', 'Nsquared']

#varnames = ['QN', 'QV']

#varnames = ['CLOUD AMOUNT']

#varnames = ['QN', 'QV']

#varnames = ['W']

#varnames = ['CLOUD AMOUNT']

#varnames = ['TABS', 'Nsquared', 'URADIAL']

#varnames = ['QV']

#varnames = ['URADIAL']

#varnames = ['THETA']

#varnames = ['N']

#varnames = ['W_ADIABATIC']

#varnames = ['DSDZ']

#varnames = ['W_DIABATIC']

#varnames = ['URADIAL']

#varnames = ['URADIAL', 'N']

varnames = ['W']

#varnames = ['Nsquared']

#MAYBE CALCULATE THETA, AND THEN STABILITY?

labels=[]
hs1=[]
hs2=[]
hs3=[]
    
for i, domsize in enumerate(domsizes):
    
    print 'domain size', domsize
    if domsize == 768:
        nave=10
    elif domsize == 1536:
        nave=10
    else:
        nave=5
        
    aveperiod2D = nave*ntave2D
    aveperiod3D = nave*ntave3D
    
    nc_inSTAT = nc_STATs[i]
    nc_in2D = nc_2Ds[i]
    nc_in3D = nc_3Ds[i]
    t2=t2s[i]
    t3=t3s[i]
    
    print 'loading 2D variables'
    nc_data2D = Dataset(nc_in2D)
    varis2D = nc_data2D.variables
    
    PW = varis2D['PW'][t2-aveperiod2D:t2,:,:]
    
    print 'loading 3D variables'
    nc_data3D = Dataset(nc_in3D)
    varis3D = nc_data3D.variables
    
    print 'loading STAT variables'
    nc_dataSTAT = Dataset(nc_inSTAT)
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
    
    
    #EDIT HIS TO CHANGE VARIABLE FOR MOIST REGION THRESHOLD (2D threshold)
    mvarname = 'PW'
    mvar = varis2D[mvarname]
    mfield = mvar[t2-aveperiod2D:t2,:,:]
    mfield_t= np.mean(mfield, axis=0)

    #EDIT THIS TO CHANGE HOW MANY STD ABOVE MEAN FOR MOIST REGION THRESHOLD
    a=1.5
    
    mfieldcrit = np.mean(mfield) + a*np.std(mfield)
    moist_points = mfield_t >= mfieldcrit
    dry_points = mfield_t < mfieldcrit 
    
    #UNCOMMENT TO USE W > W_crit AS CONVECTIVE REGION THRESHOLD
    W=varis3D['W'][t3-aveperiod3D:t3,:,:]
    W_tave = np.mean(W, axis=0)
    
        
    db=16
    W_crit = 0.01
    
    W_tave = blockave3D(W_tave, db)
    
    nxprime = nx/db
    nyprime = ny/db
    

    #moist_points = W_tave >= W_crit
    #dry_points = W_tave < W_crit
    
    #find BL to get w_BL and w_m
    dthetadz = (theta_tave[1:]-theta_tave[:-1])/np.diff(z)
    
    dthetadz_crit = 2e-3
    
    BLi = np.where(dthetadz > dthetadz_crit)[0][0]
    BLi = 0
    
    p_BL = p[BLi]
    z_BL = z[BLi]
    zBL = z[:BLi]
    p_t = 50
    z_t = z[p <= p_t*1e2][0]
    
    
    print 'z_BL (m)', z_BL

    print 'calculating dry region and convective region average profiles'
    print 'domain length = {:3.0f}, day {:3.0f} to day {:3.0f}'.format(domsize, t3D[0], t3D[-1])

    
    for k, varname in enumerate(varnames):
        print 'loading', varname
        if varname == 'URADIAL':
            #plot radial velocity in moist region/convective region?
            U_rfname = glob.glob(foutdata + '*{:d}*{:3.0f}to{:3.0f}*'.format(domsize, t3D[0], t3D[-1]))[0]
            field_tave = np.load(U_rfname)
            field_tave = blockave3D(field_tave, db)
            
        elif varname == 'BUOY':
            buoy_fname = glob.glob(foutdata + '*{:d}*_buoy_*to{:3.0f}*'.format(domsize, t3D[-1]))[0]
            print 'loading', buoy_fname
            field_tave = np.load(buoy_fname)
            field_tave = np.mean(field_tave[-nave:,:,:,:], axis=0)
            field_tave = blockave3D(field_tave, db)
        
        elif varname == 'CLOUD AMOUNT':
            QN = varis3D['QN'][t3-aveperiod3D:t3,:,:,:]
            QN_tave = np.mean(QN, axis=0)
            field_tave = blockave3D(QN_tave, db)
            
        elif varname == 'THETA':
            p0=1000*1e2
            T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
            T_tave = np.mean(T, axis=0)
            p3D = np.zeros((nz, nx, ny))
            p3D= p3D.T
            p3D[:,:,:] = p
            p3D = p3D.T
            #p3D = p3D.T
            #p4D = np.tile(p3D[:,:,:, np.newaxis], T.shape[0])
            #p4D = p4D.T
            THETA_tave = T_tave*(p0/p3D)**(c.Rd/c.cp)
            field_tave = blockave3D(THETA_tave, db)
            #T_tave = np.mean(T, axis=0)
         
        #adiabatic warming + radiative cooling balance here   
        elif varname == 'ADBWARMplusQRAD':
            
            vari = varis3D['QRAD']
            Qr = varis3D['QRAD'][t3-aveperiod3D:t3,:,:,:]
            w = varis3D['W'][t3-aveperiod3D:t3,:,:,:]
            T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
            nt = T.shape[0]
            z3D = np.zeros((ny, nx, nz))
            #z3D = z3D.T
            z3D[:,:,:] = z
            z3D = z3D[:,:,:,np.newaxis]
            z4D = np.tile(z3D,nt)
            z4D = z4D.T
            delz3D = np.zeros((ny, nx, nz-1))
            delz3D[:,:,:] = np.diff(z)
            delz3D = delz3D[:,:,:,np.newaxis]
            delz4D = np.tile(delz3D,nt)
            delz4D = delz4D.T
            gc.collect()
            
            #delz4D = np.zeros((nt, nz-1, nx, ny))
            #delz4D = delz4D.T
            #delz4D[:,:,:] = np.diff(z)
            #delz4D = delz4D.T
            s = c.cp*T + c.g*z4D
            dsdz = (s[:,1:,:,:] - s[:,:-1,:,:])/delz4D
            adbwarm = -(1./c.cp)*np.multiply(w[:,:-1,:,:], dsdz)
            adbwarm = adbwarm*(3600*24) #convert from K/s to K/day
            gc.collect()
            field_tave = adbwarm + Qr[:,:-1,:,:]
            field_tave = np.mean(field_tave, axis=0)
            field_tave = blockave3D(field_tave, db)
         
        elif varname == 'ADBWARMplusQRADBAR':
            vari = varis3D['QRAD']
            Qr = varis3D['QRAD'][t3-aveperiod3D:t3,:,:,:]
            w = varis3D['W'][t3-aveperiod3D:t3,:,:,:]
            T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
            nt = T.shape[0]
            
             #uncomment to calculate wbar*ds/dzbar (where bar is the time mean)
            w_tave = np.mean(w, axis=0)
            Qr_tave = np.mean(Qr, axis=0)
            T_tave = np.mean(T, axis=0)
             
            z3D = np.zeros((ny, nx, nz))
            z3D[:,:,:] = z
            z3D = z3D.T
            delz3D = np.zeros((ny, nx, nz-1))
            delz3D[:,:,:] = np.diff(z)
            delz3D = delz3D.T
            s = c.cp*T_tave + c.g*z3D
            dsdz = (s[1:,:,:] - s[:-1,:,:])/delz3D
            adbwarm_tave = -(1./c.cp)*np.multiply(w_tave[:-1,:,:], dsdz)
            adbwarm_tave = adbwarm_tave*(3600*24) #convert from K/s to K/day
            #adbwarm_tave = np.mean(adbwarm, axis=0) 
            field_tave = adbwarm_tave + Qr_tave[:-1,:,:]
            field_tave = blockave3D(field_tave, db)
            
            
        elif varname == 'DSDZ':
            
            T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
            T_tave = np.mean(T, axis=0)
        
            z3D = np.zeros((nz, nx, ny))
            z3D= z3D.T
            z3D[:,:,:] = z
            z3D = z3D.T
            
            s = c.cp*T_tave + c.g*z3D
            
            diffz = np.diff(z)
            diffz3D = np.zeros((nz-1, nx, ny))
            diffz3D = diffz3D.T
            diffz3D[:,:,:] = diffz
            diffz3D = diffz3D.T
            
            dsdz = (s[1:,:,:] - s[:-1,:,:])/(diffz3D)
            
            field_tave = blockave3D(dsdz, db)
            
        elif varname == 'Nsquared':
            p0=1000*1e2
            T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
            T_tave = np.mean(T, axis=0)
            p3D = np.zeros((nz, nx, ny))
            p3D= p3D.T
            p3D[:,:,:] = p
            p3D = p3D.T
            #p3D = p3D.T
            #p4D = np.tile(p3D[:,:,:, np.newaxis], T.shape[0])
            #p4D = p4D.T
            THETA_tave = T_tave*(p0/p3D)**(c.Rd/c.cp)
            diffTHETA = THETA_tave[1:,:,:] - THETA_tave[:-1,:,:]
            diffz = np.diff(z)
            diffz3D = np.zeros((nz-1, nx, ny))
            diffz3D = diffz3D.T
            diffz3D[:,:,:] = diffz
            diffz3D = diffz3D.T
            THETA_mid = (THETA_tave[1:,:,:] + THETA_tave[:-1,:,:])/2.
            N_tave = (c.g*diffTHETA)/np.multiply(diffz3D, THETA_mid)
            field_tave = blockave3D(N_tave, db)
            
        elif varname == 'W_ADIABATIC':
            #p0=1000*1e2
            T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
            T_tave = np.mean(T, axis=0)
            QRAD = varis3D['QRAD'][t3-aveperiod3D:t3,:,:,:]
            QRAD_tave = np.mean(QRAD, axis=0)
            #QRAD_mid = (QRAD_tave[1:,:,:] + QRAD_tave[:-1,:,:])/2.
        
            z3D = np.zeros((nz, nx, ny))
            z3D= z3D.T
            z3D[:,:,:] = z
            z3D = z3D.T
            
            s = c.cp*T_tave + c.g*z3D
            
            diffz = np.diff(z)
            diffz3D = np.zeros((nz-1, nx, ny))
            diffz3D = diffz3D.T
            diffz3D[:,:,:] = diffz
            diffz3D = diffz3D.T
            
            dsdz = (s[1:,:,:] - s[:-1,:,:])/(diffz3D)
            
            W_adb = (c.cp*QRAD_tave[:-1,:,:])/(dsdz*3600*24)
            field_tave = blockave3D(W_adb, db)
            
        elif varname == 'W_DIABATIC':
            T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
            T_tave = np.mean(T, axis=0)
            QRAD = varis3D['QRAD'][t3-aveperiod3D:t3,:,:,:]
            QRAD_tave = np.mean(QRAD, axis=0)
            #W = varis3D['W'][t3-aveperiod3D:t3,:,:,:]
            #W_tave = np.mean(W, axis=0)
            
            #QRAD_mid = (QRAD_tave[1:,:,:] + QRAD_tave[:-1,:,:])/2.
        
            z3D = np.zeros((nz, nx, ny))
            z3D= z3D.T
            z3D[:,:,:] = z
            z3D = z3D.T
            
            s = c.cpd*T_tave + c.g*z3D
            
            diffz = np.diff(z)
            diffz3D = np.zeros((nz-1, nx, ny))
            diffz3D = diffz3D.T
            diffz3D[:,:,:] = diffz
            diffz3D = diffz3D.T
            
            dsdz = (s[1:,:,:] - s[:-1,:,:])/(diffz3D)
            
            W_adb = (c.cp*QRAD_tave[:-1,:,:])/(dsdz*3600*24)
            
            W_diab = W_tave[:-1,:,:] - blockave3D(W_adb, db)
            field_tave = W_diab
        #calculate RH
        elif varname == 'TABS':
            varname = 'RH'
            QV = varis3D['QV'][t3-aveperiod3D:t3,:,:,:]
            QV_tave = np.mean(QV, axis=0)
            TABS = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
            TABS_tave = np.mean(TABS, axis=0)
            RH_tave = np.zeros(QV_tave.shape)
            for pi, plev in enumerate(p):
                wsat_tave = 1000*wsat(TABS_tave[pi,:,:], plev) #convert to g/kg
                RH_tave[pi,:,:] = 100*(QV_tave[pi,:,:]/wsat_tave) #convert to percent
            field_tave = RH_tave
            field_tave = blockave3D(field_tave, db)
            
        elif varname == 'VERTMASSFLUX':
            w = varis3D['W'][t3-aveperiod3D:t3,:,:,:]
            T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
            nt = T.shape[0]
            p3D = np.zeros((ny, nx, nz))
            p3D[:,:,:] = p
            p3D = p3D[:,:,:,np.newaxis]
            p4D = np.tile(p3D,nt)
            p4D = p4D.T
            gc.collect()
            rho = p4D/(T*c.Rd)
            massflux = np.multiply(w, rho)
            field_tave = np.mean(massflux, axis=0)
            field_tave = blockave3D(field_tave, db)

        else:
            vari = varis3D[varname]
            field = varis3D[varname][t3-aveperiod3D:t3,:,:,:]
            field_tave = np.mean(field, axis=0)
            field_tave = blockave3D(field_tave, db)
   

        
        #moist_field = field_tave[:,moist_points]
        #dry_field = field_tave[:,dry_points]
        
        moist_prof=[]
        dry_prof=[]
        edge_prof =[]
        
        Wedge_u = 3
        Wedge_l = 0.5
        
        print 'calculating convective/dry profile'
        for zi in np.arange(field_tave.shape[0]):
            moist_points = W_tave[zi,:,:] >= W_crit
            dry_points = W_tave[zi,:,:] < W_crit
            edge_points = np.bitwise_and(W_tave[zi,:,:] <= Wedge_u*W_crit, W_tave[zi,:,:] >= Wedge_l*W_crit)
            #print 'z', z[zi]
           # print 'len edge points', len(edge_points)
            if varname == 'CLOUD AMOUNT':
                QNz_m = field_tave[zi,moist_points]
                QNz_d = field_tave[zi,dry_points]
                QNz_edge = (field_tave[zi, edge_points])
                totpoints = 1.*nxprime*nyprime
                cldpoints_m = len(QNz_m[QNz_m > 0])
                cldpoints_d = len(QNz_d[QNz_d > 0])
                cldedgepoints = len(QNz_edge[QNz_edge > 0])
                moist_bar = (100*cldpoints_m)/totpoints
                dry_bar = (100*cldpoints_d)/totpoints
                edge_bar = (100*cldedgepoints)/totpoints
                moist_prof= moist_prof + [moist_bar]
                dry_prof = dry_prof + [dry_bar]
                edge_prof = edge_prof + [edge_bar]
            else:    
                moist_bar = np.mean(field_tave[zi, moist_points])
                dry_bar = np.mean(field_tave[zi, dry_points])
                edge_bar = np.mean(field_tave[zi, edge_points])
                #print edge_bar
                moist_prof= moist_prof + [moist_bar]
                dry_prof = dry_prof + [dry_bar]
                edge_prof = edge_prof + [edge_bar]
                
        moist_prof = np.array(moist_prof)
        dry_prof = np.array(dry_prof)
        edge_prof = np.array(edge_prof)
            
        #moist_prof = np.mean(moist_field, axis=1)
        #dry_prof = np.mean(dry_field, axis=1)
    
        if varname == 'RH':
            units='%'
        elif varname == 'URADIAL':
            units='m/s'
        elif varname == 'THETA':
            units='K'
        elif varname == 'CLOUD AMOUNT':
            units = '%'
        elif varname == 'Nsquared':
            units = r's$^{{-1}}$'
        elif varname == 'W_ADIABATIC':
            units = 'm/s'
        elif varname == 'W_DIABATIC':
            units = 'm/s'
        elif varname == 'DSDZ':
            units = 'J/m'
        elif varname == 'VERTMASSFLUX':
            units = r'kg m$^{{-2}}$ s$^{{-1}}$'
            #cb.set_label(units)
        elif varname == 'BUOY':
            units = 'Pa'
        else:
            units = vari.units.strip()
            
        if varname == 'ADBWARMplusQRAD':
            titlename = r'$-w\frac{ds}{dz}$ + Q$_r$'
        elif varname == 'ADBWARMplusQRADBAR':
            titlename =  r'$-\bar{w}\frac{d\bar{s}}{dz}$ + Q$_r$'
        elif varname == 'QRAD':
            titlename = r'Q$_r$'
        elif varname == 'ADBWARM':
            titlename = r'$\-wfrac{ds}{dz}$'
        elif varname == 'TABSprime':
            titlename = r'$T - \bar{T}$'
        elif varname == 'DSDZ':
            titlename = r'$\frac{ds}{dz}$'
        elif varname == 'Nsquared':
            titlename = r'$N^2$'
        elif varname == 'VERTMASSFLUX':
            titlename = r'$\rho w$'
        elif varname == 'URADIAL':
            titlename = r'$u_r$'
        elif varname == 'BUOY':
            titlename = r'$\rho \frac{dw}{dt}$'
        else:
            titlename = varname
            
        z_temp = z[:moist_prof.shape[0]]
            
        interior = np.bitwise_and(z_temp > z_BL, z_temp < z_t)
        zint = z_temp[interior]
        moist_prof = moist_prof[interior] 
        dry_prof = dry_prof[interior]
        edge_prof = edge_prof[interior]
        
        moist_profave = np.mean(moist_prof[np.isfinite(moist_prof)])
        dry_profave = np.mean(dry_prof[np.isfinite(dry_prof)])
        #edge_profave = np.mean(edge_prof)
        
        print 'moist region vertical average {:s} = {:4.4f} ({:s})'.format(varname, moist_profave, units)
        print 'dry region vertical average {:s} = {:4.4f} ({:s})'.format(varname, dry_profave, units)
        
        fig = plt.figure(k)
        ax1 = fig.add_subplot(2,1,1)
        #plt.suptitle(r'{:s} ({:s}), convective region {:s} threshold = {:2.1f} {:s}'.format(varname, units, mvarname, mfieldcrit, mvar.units.strip()))
        h1, = ax1.plot(moist_prof, zint/1e3, color=colors[i])
        if db == 1:
            ax1.set_title('{:s}, convective region using W > W$_c$ = {:3.3f} m/s'.format(titlename, W_crit))
        else:
            ax1.set_title('{:s}, convective region using W > W$_c$ = {:3.3f} m/s, block-averaging over ({:2.0f} km)$^2$'.format(titlename, W_crit, db*(np.diff(x)[0])/1e3))
        #axarr[0,].set_xlabel('{:s} ({:s})'.format(varname, vari.units.strip()))
        ax1.set_ylabel('z (km)')
        ax1.set_ylim(z_BL/1e3, z_t/1e3)
        #if (varname != 'QV'):
        #    axarr[0,].axvline(moist_profave, color='k', alpha=0.5, label='vertical average')
        ax2 = fig.add_subplot(2,1,2)
        ax2.plot(dry_prof, zint/1e3, color=colors[i])
        if db == 1:
            ax2.set_title('{:s}, convection-free region using W < W$_c$ = {:3.3f} m/s'.format(titlename, W_crit))
        else:
            ax2.set_title('{:s}, convection-free region using W < W$_c$ = {:3.3f} m/s, block-averaging over ({:2.0f} km)$^2$'.format(titlename, W_crit, db*(np.diff(x)[0])/1e3))
        ax2.set_xlabel('{:s} ({:s})'.format(titlename, units))
        ax2.set_ylabel('z (km)')
        ax2.set_ylim(z_BL/1e3, z_t/1e3)
        if varname == 'W_ADIABATIC':
            ax1.set_xlim((-0.05, 0.16))
            ax2.set_xlim((-0.008, 0))
        if varname == 'W':
            ax1.set_xlim(0, 0.2)
            #ax1.set_xlim(0.5, 0.8)
        if varname == 'RH':
            ax1.set_xlim(0, 100)
            ax2.set_xlim(0, 100)
        #if varname == 'W_DIABATIC':
        #   ax1.set_xlim((
        #axarr[1,].axvline(dry_profave, color='k', alpha=0.5, label='vertical average')
        label1=r'day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$'.format(t3D[0], t3D[-1], domsize)
        #fig.legend(hs, labels, 'best')       
        fig.savefig(fout + '{:s}dryconvprof_Wthresh_blkavg_day250_{:d}day_db{:d}.pdf'.format(varname, nave, db))
        fig = plt.figure(len(varnames)+k)
        ax = fig.gca()
        h2, = ax.plot(moist_prof - dry_prof, zint/1e3, color=colors[i])
        ax.set_title('differential {:s} between convective region and dry region'.format(varname, nave))
        ax.set_ylabel('z (km)')
        ax.set_xlabel('differential {:s} ({:s})'.format(varname, units))
        ax.set_ylim(z_BL/1e3, z_t/1e3)
        fig.savefig(fout + 'DIFF{:s}dryconvprof_Wthresh_blkavg_day250_{:d}day_db{:d}.pdf'.format(varname, nave, db))
        
        if varname == 'URADIAL':
            fig = plt.figure(2*len(varnames)+k)
            ax = fig.gca()
            h3, = ax.plot(edge_prof, zint/1e3, color=colors[i])
            ax.set_title(r'{:s} at convective edge, {:3.1f}W$_c$  < W < {:3.1f}W$_c$  m/s'.format(varname, Wedge_l, Wedge_u))
            ax.set_ylabel('z (km)')
            ax.set_xlabel('{:s} ({:s})'.format(varname, units))
            ax.set_ylim(z_BL/1e3, z_t/1e3)
            fig.savefig(fout + '{:s}convedge_Wthresh_blkavg_day250_{:d}day_db{:d}.pdf'.format(varname, nave, db))
            
    labels = labels + [label1]
    hs1 = hs1 + [h1]
    hs2 = hs2 + [h2] 
    if varname == 'URADIAL':
        hs3 = hs3 + [h3]
        
    gc.collect()

#set up legends
for k, varname in enumerate(varnames):
    fig=plt.figure(k)
    ax=plt.subplot(2,1,1)
    ax.legend(hs1, labels, loc='best')
    fig.savefig(fout + '{:s}dryconvprof_Wthresh_blkavg_day250_{:d}day_db{:d}.pdf'.format(varname, nave, db))
    #plt.close(k)
    fig=plt.figure(len(varnames)+k)
    ax = fig.gca()
    ax.legend(hs2, labels, loc='best')
    fig.savefig(fout + 'DIFF{:s}dryconvprof_Wthresh_blkavg_day250_{:d}day_db{:d}.pdf'.format(varname, nave, db))
    plt.close(len(varname)+k)
    if varname == 'URADIAL':
        fig=plt.figure(2*len(varnames)+ k)
        ax = fig.gca()
        ax.legend(hs3, labels, loc='best')
        fig.savefig(fout + '{:s}convedge_Wthresh_blkavg_day250_{:d}day_db{:d}.pdf'.format(varname, nave, db))
        plt.close(2*len(varnames)+k)
        


    
        
    ##simple model results   
    #w_BL = -.0066
    #w_m = 0.0334
    #
    #W_m = 3
    #W_d = -4.5
    #
    #q_FA = 0.000175*1e3 #in g/kg
    #
    #T_s=302
    #h=7
    #gamma_PH=6.2
    #zeta_T=16
    #p_s = 1000e2
    #
    #def findT(T_s, p):
    #    
    #    zeta = -h*(np.log(p/p_s))
    #    
    #    if (zeta < zeta_T):
    #        T = T_s - gamma_PH*zeta
    #    else:
    #        T = T_s - gamma_PH*zeta_T
    #    return T
    #    
    #qv_s = np.zeros(p.size)
    #    
    #for i, plev in enumerate(p):
    #    qv = wsat(findT(T_s, plev), p_s)
    #    qv_s[i] = qv*1e3    
    #    
    #qv_s = qv_s[interior]

    
    #if varname == 'W':
    #    axarr[1,].axvline(w_BL, color='r', alpha=0.5, label='simple model')
    #    axarr[0,].axvline(w_m, color='r', alpha=0.5, label='simple model')
    #if varname == 'QRAD':
    #    axarr[0,].axvline(W_m, color='r', label = 'simple model')
    #    axarr[1,].axvline(W_d, color='r', label = 'simple model')
    #    axarr[0,].set_xlim(-4, 3.5)
    #if varname == 'QV':
    #    axarr[1,].axvline(q_FA, color='r', alpha=0.5, label = 'simple model')
    #    axarr[0,].plot(qv_s, color='r', alpha=0.5, label = 'simple model')
    #
    

















