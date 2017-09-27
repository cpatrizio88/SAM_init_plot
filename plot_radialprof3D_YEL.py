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
from misc_fns import radprof3D, radprof
#from SAM_init_plot.block_fns import blockave2D, blockave3D, blockxysort2D
from block_fns import blockave2D, blockave3D, blockxysor2D
#from thermolib.constants import constants
from constants import constants
from wsat import wsat
#from thermolib.wsat import wsat
import gc

c=constants()

matplotlib.rcParams.update({'font.size': 26})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 2})
matplotlib.rcParams.update({'legend.fontsize': 24})
matplotlib.rcParams.update({'mathtext.fontset': 'cm'})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fdata = '/Users/cpatrizio/data/SST302/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_RADIALXSECTIONNEW/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to195_64vert_ubarzero_RADIALXSECTIONNEW/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to190_64vert_ubarzero_RADIALXSECTIONNEW/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_RADIAL/'

fpath2D = '/glade/scratch/patrizio/OUT_2D_nc/'
fpath3D = '/glade/scratch/patrizio/OUT_3D_nc/'
fout = '/glade/scratch/patrizio/OUT_3D_FIGS/'

#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*day230to250*302K.nc')[0]

#nc_in2D = glob.glob(fpath2D + '*512x512*3000m*day180to195*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day180to195*302K.nc')[0]

#nc_in2D = glob.glob(fpath2D + '*1024x1024*3000m*day170to180*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]

nc_in2D = glob.glob(fpath2D + '*2048*3000m*.nc')[0]
nc_in3D = glob.glob(fpath3D + '*2048*3000m*.nc')[0]

#domsize=768
#domsize=1536
domsize=3072

nave=10
ntave2D=24
ntave3D=4
t2 = -35*ntave2D
t3 = -35*ntave3D
t2= -1
t3= -1

print 'domsize', domsize

aveperiod2D = nave*ntave2D
aveperiod3D = nave*ntave3D

print 'loading .nc files'
nc_data3D = Dataset(nc_in3D)
nc_data2D = Dataset(nc_in2D)
varis3D = nc_data3D.variables
varis2D = nc_data2D.variables

print 'loading 2D variables'
nc_data2D = Dataset(nc_in2D)
varis2D = nc_data2D.variables

PW = varis2D['PW'][t2-aveperiod2D:t2,:,:]

print 'loading 3D variables'
nc_data3D = Dataset(nc_in3D)
varis3D = nc_data3D.variables

t3D = varis3D['time'][t3-aveperiod3D:t3]
t2D = varis2D['time'][t2-aveperiod2D:t2]
x = varis3D['x'][:]
y = varis3D['y'][:]
z = varis3D['z'][:]
p = varis3D['p'][:]
p = p*1e2

#varnames = ['QRAD', 'W', 'U', 'QN', 'QV']
#varnames=['QV', 'TABS']
#varnames=['TABS', 'QN', 'QRAD', 'W']
#varnames = ['W' ]
varnames = ['QRAD']
#varnames = ['W']
#varnames=['dWdt']
#varnames=['BuoyFlux']
#varnames=['Nsquared']
#varnames = ['DSDZ']

varnames = ['DSDZ', 'Nsquared', 'QN', 'QP', 'RELH']
varnames = ['TABS', 'W']
varnames = ['VERTMASSFLUX']
varnames = ['W']
varnames = ['URADIAL']
varnames = ['U']

#varnames = ['QN', 'TABS', 'W']
#varnames = ['QP']
varnames = ['QN', 'QP', 'TABS']
#varnames = ['QP']
#varnames = ['QRAD', 'W', 'TABS']
#varnames = ['DSDZ', 'Nsquared', 'URADIAL', 'W']
#varnames = ['URADIAL', 'QRAD']
#varnames = ['VERTMASSFLUX', 'URADIAL']
#varnames = ['BUOY', 'BuoyFlux']
varnames = ['BUOY', 'Nsquared', 'TABS']
varnames = ['TABS']
varnames = ['BuoyFlux']
varnames = ['W']
varnames = ['TABS']
varnames = ['BuoyFlux']

varnames = ['W']
#varnames = ['TABS', 'Nsquared']
varnames = ['TABS']


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
db=1
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

zedges = np.zeros(z.shape)
zedges[1:] = (z[1:] + z[:-1])/2.

nbins=[rbins, zedges]

print 'calculating z-r plots'
for varname in varnames:
    print 'loading', varname
    if varname == 'URADIAL':
        fname = glob.glob(fdata + '{:d}*URADIAL_*to{:3.0f}*'.format(domsize, t3D[-1]))[0]
        print 'loading', fname
        field_tave = np.load(fname)
    elif varname == 'dWdt':
        fname = glob.glob(fdata + '{:d}*dwdt*'.format(domsize))[0]
        print 'loading', fname
        dwdt = np.load(fname)
        #dwdt already daily averaged
        dwdt = dwdt[t3/ntave3D,:,:,:]
        field_tave = dwdt
    elif varname == 'BuoyFlux':
        fname = glob.glob(fdata + '{:d}*buoyflux_*to{:3.0f}*'.format(domsize, t3D[-1]))[0]
        print 'loading', fname
        buoyflux = np.load(fname)
        #dwdt already daily averaged
        buoyflux = np.mean(buoyflux[t3/ntave3D-nave:t3/ntave3D,:,:,:],axis=0)
        field_tave = buoyflux
        delz3D = np.zeros((nx, ny, nz-1))
        delz3D[:,:,:] = np.diff(z)
        delz3D = delz3D.T
        field_tave = field_tave*delz3D
        znew = z[:-1]
    elif varname == 'BUOY':
        buoy_fname = glob.glob(fdata + '*{:d}*_buoy_*to{:3.0f}*'.format(domsize, t3D[-1]))[0]
        print 'loading', buoy_fname
        field_tave = np.load(buoy_fname)
        field_tave = np.mean(field_tave[t3/ntave3D-nave:t3/ntave3D,:,:,:],axis=0)
        delz3D = np.zeros((nx, ny, nz-1))
        delz3D[:,:,:] = np.diff(z)
        delz3D = delz3D.T
        field_tave = field_tave*delz3D
        znew = z[:-1]
    #calculate relative humidity 
    elif varname == 'TABS':
        varname = 'RH'
        QV = varis3D['QV'][t3-aveperiod3D:t3,:,:,:]
        T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
        T_tave = np.mean(T, axis=0)
        QV_tave = np.mean(QV, axis=0)
        RH_tave = np.zeros(QV_tave.shape)
        for i, plev in enumerate(p):
            wsat_tave = 1000*wsat(T_tave[i,:,:], plev) #convert to g/kg
            RH_tave[i,:,:] = 100*(QV_tave[i,:,:]/wsat_tave) #convert to percent
        field_tave = RH_tave
        
        ##UNCOMMENT TO CALCULATE TABS PERTURBATION
        #vari = varis3D['TABS']
        #varname = 'TABSprime'
        #TABS = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
        #T_tave = np.mean(TABS, axis=0)
        #T_bar = np.mean(np.mean(T_tave, axis=2), axis=1)
        #T_bar3D = np.zeros(T_tave.shape)
        #T_bar3D = T_bar3D.T
        #T_bar3D[:,:,:] = T_bar
        #T_bar3D = T_bar3D.T
        #Tprime = T_tave - T_bar3D
        #field_tave = Tprime
    #calculate wind speed
    elif varname == 'U':
       vari = varis3D['U']
       U = varis3D['U'][t3-aveperiod3D:t3,:,:,:]
       V = varis3D['V'][t3-aveperiod3D:t3,:,:,:]
       V_tave = np.mean(V, axis=0)
       U_tave = np.mean(U, axis=0)
       #V = V[ntrunc:,:,:,:]
       #U = varis3D['U'][t3-aveperiod3D:t3,:,:,:]
       #U = U[ntrunc:,:,:,:]
       field_tave = np.sqrt(np.power(V_tave, 2) + np.power(U_tave, 2))
    #field_tmp = field.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
    #field_tave = np.mean(field_tmp, axis=1)
    #calculate adiabatic warming + QRAD balance 
    #elif varname == 'QRAD':
    #   # #uncomment to calculate wbar*ds/dzbar (where bar is the time mean)
    #   #gc.collect()
    #   #vari = varis3D['QRAD']
    #   #varname = 'ADBWARMplusQRADBAR'
    #   #Qr = varis3D['QRAD'][t3-aveperiod3D:t3,:,:,:]
    #   #Qr_tave = np.mean(Qr, axis=0)
    #   #w = varis3D['W'][t3-aveperiod3D:t3,:,:,:]
    #   #T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
    #   #
    #   #w_tave = np.mean(w, axis=0)
    #   #
    #   #T_tave = np.mean(T, axis=0)
    #   # 
    #   #z3D = np.zeros((ny, nx, nz))
    #   #z3D[:,:,:] = z
    #   #z3D = z3D.T
    #   #delz3D = np.zeros((ny, nx, nz-1))
    #   #delz3D[:,:,:] = np.diff(z)
    #   #delz3D = delz3D.T
    #   #s = c.cp*T_tave + c.g*z3D
    #   #dsdz = (s[1:,:,:] - s[:-1,:,:])/delz3D
    #   #adbwarm_tave = -(1./c.cp)*np.multiply(w_tave[:-1,:,:], dsdz)
    #   #adbwarm_tave = adbwarm_tave*(3600*24) #convert from K/s to K/day
    #   ##adbwarm_tave = np.mean(adbwarm, axis=0)
    #   #field_tave = adbwarm_tave + Qr_tave[:-1,:,:]
    #   
    #   
    #   
    #   gc.collect()
    #   vari = varis3D['QRAD']
    #   #varname = 'ADBWARMplusQRADBAR'
    #   varname = 'ADBWARMplusQRAD'
    #   Qr = varis3D['QRAD'][t3-aveperiod3D:t3,:,:,:]
    #   #Qr_tave = np.mean(Qr, axis=0)
    #   w = varis3D['W'][t3-aveperiod3D:t3,:,:,:]
    #   T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
    #   nt = T.shape[0]
    #   z3D = np.zeros((ny, nx, nz))
    #   #z3D = z3D.T
    #   z3D[:,:,:] = z
    #   z3D = z3D[:,:,:,np.newaxis]
    #   z4D = np.tile(z3D,nt)
    #   z4D = z4D.T
    #   delz3D = np.zeros((ny, nx, nz-1))
    #   delz3D[:,:,:] = np.diff(z)
    #   delz3D = delz3D[:,:,:,np.newaxis]
    #   delz4D = np.tile(delz3D,nt)
    #   delz4D = delz4D.T
    #   gc.collect()
       
       #delz4D = np.zeros((nt, nz-1, nx, ny))
       #delz4D = delz4D.T
       #delz4D[:,:,:] = np.diff(z)
       #delz4D = delz4D.T
       s = c.cp*T + c.g*z4D
       dsdz = (s[:,1:,:,:] - s[:,:-1,:,:])/delz4D
       adbwarm = -(1./c.cp)*np.multiply(w[:,:-1,:,:], dsdz)
       adbwarm = adbwarm*(3600*24) #convert from K/s to K/day
       field_temp = adbwarm + Qr[:,:-1,:,:] 
       field_tave = np.mean(field_temp, axis=0)
       znew = z[:-1]
    #elif varname == 'W':
    #   
    #   varname = 'ADBWARM'
    #   w = varis3D['W'][t3-aveperiod3D:t3,:,:,:]
    #   T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
    #   #w_t = np.mean(w[t-nave3D:t,:,:,:], axis=0)
    #   w_tave = np.mean(w, axis=0)
    #   #T_t = np.mean(T[t-nave3D:t,:,:,:], axis=0)
    #   T_tave = np.mean(T, axis=0)
    #   #z3D = np.zeros(field_t.shape)
    #   z3D = np.zeros(T_tave.shape)
    #   z3D = z3D.T
    #   z3D[:,:,:] = z
    #   z3D = z3D.T
    #   delz3D = np.zeros((nz-1, nx, ny))
    #   delz3D = delz3D.T
    #   delz3D[:,:,:] = np.diff(z)
    #   delz3D = delz3D.T
    #   s = c.cp*T_tave + c.g*z3D
    #   dsdz = (s[1:,:,:] - s[:-1,:,:])/delz3D
    #   adbwarm = -(1./c.cp)*np.multiply(w_tave[:-1], dsdz)
    #   adbwarm = adbwarm*(3600*24) #convert from K/s to K/day
    ##   field_tave = adbwarm
       
    #ADD MASS FLUX (KG per S) HERE
    
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
       
    
    elif varname == 'DSDZ':
       T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
       T_tave = np.mean(T, axis=0)
       #z3D = np.zeros(field_t.shape)
       z3D = np.zeros(T_tave.shape)
       z3D = z3D.T
       z3D[:,:,:] = z
       z3D = z3D.T
       delz3D = np.zeros((nz-1, nx, ny))
       delz3D = delz3D.T
       delz3D[:,:,:] = np.diff(z)
       delz3D = delz3D.T
       s = c.cp*T_tave + c.g*z3D
       dsdz = (s[1:,:,:] - s[:-1,:,:])/delz3D
       field_tave = dsdz
        
    elif varname == 'Nsquared':
        p0=1000*1e2
        T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
        T_tave = np.mean(T, axis=0)
        p3D = np.zeros((T_tave.shape))
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
        field_tave = N_tave
    else:
        vari = varis3D[varname]
        field = varis3D[varname][t3-aveperiod3D:t3,:,:,:]
        field_tave = np.mean(field, axis=0)
       
    nzplot = field_tave.shape[0]
    nxplot = field_tave.shape[1]
    nyplot = field_tave.shape[2]
    
       
    #fieldmeans has shape (redges.size, zedges.size) = nbins
    print '2d contouring'
    redges, zedges, fieldmeans = radprof3D(field_tave, xx, yy, z[:nzplot], mcenter, nbins=nbins)
    rbin_centers = (redges[1:] + redges[:-1])/2.
    zbin_centers = (zedges[1:] + zedges[:-1])/2.
    
    #rr has shape (zedges.size, redges.size) = nbins.T
    rr, zz = np.meshgrid(rbin_centers, zbin_centers)
    
    plt.figure(2)
    ax=plt.gcf().gca()
    vmin=np.ma.masked_invalid(fieldmeans).min()
    vmax=np.ma.masked_invalid(fieldmeans).max()
    #EDIT: COLORBAR FORMATTING. IF NOT MONOTONIC, SWITCH TO DIVERGING COLORBAR, OTHERWISE
    #USE SEQUENTIAL COLORBAR. SET COLOR MAPPING LIMITS MANUALLY. 
    if (vmin < 0 and vmax > 0):
       if np.abs(vmin) < np.abs(vmax):
          vmin = vmin - (vmax - np.abs(vmin))
       else:
          vmax = vmax + (np.abs(vmin) - vmax)
       cmap = cm.RdBu_r
    else:
       cmap = cm.YlGnBu
    if varname == 'W':
       #vmax=0.3
       #vmin=0
       vmin=-.01
       vmax=0
       #cmap = cm.RdBu_r #use cm.copper for subsidence?
       cmap = cm.RdYlBu_r
       #cmap = cm.afmhot
    if varname == 'VERTMASSFLUX':
       #vmin=0
       #vmax=0.15
       vmin = -.01
       vmax = 0
       cmap=cm.RdYlBu_r
    if varname == 'BUOY':
       vmin = -12
       vmax = 12
    if (varname == 'ADBWARMplusQRAD') or (varname == 'ADBWARMplusQRADBAR'):
       vmin=-6
       vmax=6
       cmap=cm.RdBu_r
    if varname == 'QRAD':
       vmin=-3
       vmax=3
       cmap=cm.RdBu_r
    if varname == 'U':
       vmin=0
       vmax=9
       cmap = cm.RdYlBu_r
    if varname == 'QN':
       vmin=0
       vmax=0.35
       cmap = cm.RdYlGn_r
    if varname == 'QP':
       #vmin=0
       #vmax=0.35
       vmin=0
       vmax=0.01
       cmap = cm.RdYlGn_r
    if varname == 'Nsquared':
       vmin = 1.5e-4
       vmax = 2.7e-4
       cmap = cm.YlGnBu
    if varname == 'DSDZ':
       vmin = 0
       vmax = 15
       cmap = cm.YlGnBu
    if varname == 'TABSprime':
       vmin=-2
       vmax = 2
       cmap = cm.RdBu_r
    if varname == 'URADIAL':
       vmin=-6
       vmax=6
    if varname == 'BuoyFlux':
       vmin=-0.05
       vmax=0.05
       #vmin=-1
       #vmax=1
    if varname == 'RH':
       vmin=0
       vmax=100
       cmap = cm.RdYlGn_r
       
    if varname == 'ADBWARMplusQRAD':
        titlename = r'$-w\frac{ds}{dz}$ + Q$_r$'
    elif varname == 'ADBWARMplusQRADBAR':
        titlename =  r'$-[w][\frac{ds}{dz}]$ + Q$_r$'
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
        #titlename = r'$\rho \Delta z \frac{dw}{dt}$' 
        titlename = r'$\rho \Delta z \theta^\prime_v$'
    elif varname == 'U':
        titlename = r'$(u^2 + v^2)^{\frac{1}{2}}$'
    elif varname == 'BuoyFlux':
        #titlename = r'$\rho \Delta z \frac{dw^{2}}{dt}$' 
        titlename = r'$\rho w\Delta z\theta^\prime_v$'
    elif varname == 'W':
        titlename = r'$w$'
    elif varname == 'QV':
        titlename = r'$q_v$'
    elif varname == 'QN':
        titlename = '$q_n$'
    elif varname == 'QP':
        titlename = r'$q_p$'
    else:
        titlename = varname
       
    print 'plotting'
    fieldmeans = np.transpose(fieldmeans)
    #extent=[redges[0]/(domsize*1e3), redges[-1]/(domsize*1e3), zedges[0]/1e3, zedges[-1]/1e3]
    #plt.contourf(rr/(1e3*domsize), zz/(1e3), np.transpose(fieldmeans), cmap=cm.RdYlBu_r)
    #plt.imshow(np.transpose(fieldmeans), origin='lower', interpolation='nearest', aspect='auto', extent=extent, cmap=cm.RdYlBu_r) 
    if varname == 'QV':
         plt.pcolormesh(rr/(1e3*domsize), zz/(1e3), fieldmeans, vmin=vmin, vmax=vmax, cmap=cmap, norm=matplotlib.colors.LogNorm())
    else:
         plt.pcolormesh(rr/(1e3*domsize), zz/(1e3), fieldmeans, vmin=vmin, vmax=vmax, cmap=cmap)
    #plt.xlabel('radial distance (km)')
    plt.xlabel(r'$\hat{r}$', fontsize = 38)
    #plt.xlabel('fractional distance from moist region center, relative to domain size')
    plt.ylabel('z (km)', fontsize=32)
    cb=plt.colorbar()
    if varname == 'RH':
        units='%'
        #cb.set_label(units)
    elif varname == 'dWdt':
        units=r'm/s$^2$'
        #cb.set_label(units)
    elif varname == 'BuoyFlux':
        units=r'W/m$^2$'
        #cb.set_label(units)
    elif varname == 'Nsquared':
        units=r's$^{-2}$'
        #cb.set_label(units)
    elif varname == 'DSDZ':
        units = 'J/m'
        #cb.set_label(units)
    elif varname == 'VERTMASSFLUX':
        units = r'kg m$^{{-2}}$ s$^{{-1}}$'
        #cb.set_label(units)
    elif varname == 'BUOY':
        units = 'Pa'
    elif varname == 'URADIAL':
        units = 'm/s'
    else:
        units = vari.units.strip()
        #cb.set_label(units)
    cb.set_label(units)
    ax.set_ylim(0,16)
    ax.set_xlim([0, (1/np.sqrt(2))])
    #plt.suptitle('{:s}, day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$'.format(titlename, t3D[0], t3D[-1], domsize))
    tt1 = plt.title('{:s}, day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$'.format(titlename, t3D[0], t3D[-1], domsize))
    tt1 = plt.title('{:s}, day {:3.0f} to {:3.0f} average'.format(titlename, t3D[0], t3D[-1]))
    tt1 = plt.title('{:s}, domain size = ({:d} km)$^2$'.format(titlename, domsize))
    if (varname == 'ADBWARMplusQRAD') or (varname == 'ADBWARMplusQRADBAR') or (varname == 'BUOY') or (varname == 'BuoyFlux'):
            tt1.set_position([0.5, 1.04])
    else:
        tt1.set_position([0.5, 1.02])
    #plt.title('{:s} ({:s}), day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$, x bin width = {:2.2f}'.format(varname, vari.units.strip(), t3D[t-nave*ntave3D], t3D[t], domsize, 1./nbins[0]))
    plt.savefig(fout + '{:s}radialxsectiondist_day{:3.0f}to{:3.0f}_{:d}dayNEW.pdf'.format(varname, t3D[0], t3D[-1], nave))
    #plt.savefig(fout + '{:s}radialxsection_day{:3.0f}to{:3.0f}.pdf'.format(varname, t3D[t-nave*ntave3D], t3D[t]))
    plt.close()
    
    
    
        
    #if varname == 'W':
    ##z_BLi = np.where(z > 1.5e3)[0][0]
    #z_BLi=0
    #z_ti = np.where( z > 15e3)[0][0]
    #
    #field_tave_zave = np.mean(field_tave[z_BLi:z_ti], axis=0)
    #redges_zave, fieldmeans_zave = radprof(field_tave_zave, xx, yy, mcenter, nbins=nbins[0])
    #rcenters_zave = (redges_zave[1:] + redges_zave[:-1])/2
    #
    #plt.figure(3)
    #plt.plot(rcenters_zave/(1e3*domsize), fieldmeans_zave, 'gx', mew=2, label='{:d} km, day {:3.0f} to {:3.0f} mean'.format(domsize, t3D[0], t3D[-1]))
    ##ax.set_xlim([0, (1/np.sqrt(2))])
    #plt.xlabel(r'$\hat{r}$')
    #plt.ylabel('{:s} ({:s})'.format(titlename, units))
    #if varname == 'W':
    #    plt.ylim((-.006, 0))
    #if (varname == 'ADBWARMplusQRAD') or (varname == 'ADBWARMplusQRADBAR'):
    #    plt.ylim(-2, 0.5)
    #if varname == 'VERTMASSFLUX':
    #    plt.ylim((-.005, 0))
    #    #plt.ylim((0, 0.15))
    ##plt.ylim((-.01, 0.2))
    ##plt.axhline(0, color='k', alpha=0.4)
    #plt.title('tropospheric average {:s}'.format(titlename))
    ##plt.title('vertically-averaged {:s} ({:s}), day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$, x bin width = {:2.2f} km'.format(varname, units, t3D[0], t3D[-1], domsize, binwidth/1e3))
    #plt.savefig(fout + '{:s}radialzavedist_day{:3.0f}to{:3.0f}_{:d}dayNEW.pdf'.format(varname, t3D[0], t3D[-1], nave))
    ##plt.close()


