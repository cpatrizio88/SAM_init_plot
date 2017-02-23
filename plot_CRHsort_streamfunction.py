from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from thermolib.wsat import wsat
from thermolib.constants import constants
import SAM_init_plot.block_fns
from SAM_init_plot.block_fns import blocksort3D

c = constants()

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to130_64vert_ubarzero_CRHSORTSTREAMFN/'

#nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day90to130*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*512x512*3000m*day90to130*302K.nc')[0]

nc_in3D = glob.glob(fpath3D + '*1024x1024*3000m*day90to130*302K.nc')[0]
nc_in2D = glob.glob(fpath2D + '*1024x1024*3000m*day90to130*302K.nc')[0]
 
domsize=1536

#fraction of a day to average over.
fac=1
ndays=1

ntave2D=int(24*fac)
ntave3D=int(4*fac)

#EDIT: Time in days to look at streamfunction (t2 and t3 should correspond to the same time)
ts = np.arange(-38/fac, 0)

ts = ts.astype(int)

#ts=[50]

for t in ts:

    t2 = t*ntave2D
    #t2 = -1
    t3 = t*ntave3D
    #t3 = -1
    
    aveperiod2D = ndays*ntave2D
    aveperiod3D = ndays*ntave3D
    
    print 'loading netCDF files'
    nc_data3D = Dataset(nc_in3D)
    nc_data2D = Dataset(nc_in2D)
    varis3D = nc_data3D.variables
    varis2D = nc_data2D.variables
    
    print 'loading variables'
    t3D = varis3D['time'][t3-aveperiod3D:t3]
    t2D = varis2D['time'][t2-aveperiod2D:t2]
    print 't2D', t2D
    print 't3D', t3D
    x = varis3D['x'][:]
    y = varis3D['y'][:]
    z = varis3D['z'][:]
    p = varis3D['p'][:]
    p = p*1e2
    
    nt3D = t3D.size
    nt2D = t2D.size
    nz = z.size
    nx = x.size
    ny = y.size
    
    #calculate CRH
    PW = varis2D['PW'][t2-aveperiod2D:t2,:,:]
    SWVP = varis2D['SWVP'][t2-aveperiod2D:t2,:,:]
    CRH = PW/SWVP
    
    #need w to calculate streamfuncion
    w = varis3D['W'][t3-aveperiod3D:t3,:,:,:]
    
    #need T to calculate density
    T = varis3D['TABS'][t3-aveperiod3D:t3,:,:,:]
    
    rhoz = p/(c.Rd*np.mean(np.mean(T, axis=3), axis=2))
    
    print 'calculating daily averages'
    #reshape to do daily averages
    #T_tmp = T.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
    #w_tmp = w.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
    #rho_tmp = rhoz.reshape(nt3D/ntave3D, ntave3D, nz)
    #CRH_tmp = CRH.reshape(nt2D/ntave2D, ntave2D, nx, ny)
    
    #daily average
    #CRH_tave = np.mean(CRH_tmp, axis=1)
    #T_tave = np.mean(T_tmp, axis=1)
    #w_tave = np.mean(w_tmp, axis=1)
    #rhoz_tave = np.mean(rho_tmp, axis=1)
    
    CRH_tave = np.mean(CRH, axis=0)
    T_tave = np.mean(T, axis=0)
    w_tave = np.mean(w, axis=0)
    rhoz_tave = np.mean(rhoz, axis=0)
    
    #times in days
    #times = np.arange(int(t3D[0]), int(t3D[-1]))
    
    
    #number of grid points in a block
    db=16
    
    ###### USER EDIT: time to look at streamfunction
    #t = 0
    
    #sort CRH and w
    CRH_wsort = blocksort3D(CRH_tave[:,:], w_tave[:,:,:], db)
    CRHsort = CRH_wsort[0]
    wsort = CRH_wsort[1]
    
    #calculate the streamfunction
    print 'calculating streamfunction'
    psi = np.zeros((z.size, CRHsort.size))
    for i in range(z.size):
        for j in range(1, CRHsort.size):
            psi[i,j] = psi[i, j-1] + rhoz_tave[i]*wsort[i,j]
        
    CRHranks = np.arange(np.size(CRHsort))
    CRHrankss, zz = np.meshgrid(CRHranks, z)
    
    ########USER EDIT: variable to overlay with streamfunction
    varname = 'QN'
    vari = varis3D[varname]
    field = vari[t3-aveperiod3D:t3,:,:,:]
    #var_tmp = vari[:].reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
    field_tave = np.mean(field, axis=0)
    
    CRH_varsort = blocksort3D(CRH_tave[:,:], field_tave[:,:,:], db)
    varsort = CRH_varsort[1]
    
    if varname=='QV':
       cvals = np.arange(0, 19)
    elif varname=='QN':
       cvals = np.linspace(0, 0.2, 20)
    else:
       cvals=20
    
    plt.figure()
    ax=plt.gcf().gca()
    plt.contour(CRHrankss, zz/1e3, psi, 15, linewidths=1.5)
    plt.contourf(CRHrankss, zz/1e3, varsort, cvals, cmap=cm.RdYlBu_r)
    ax.set_ylim((0,15))
    plt.xlabel('CRH rank')
    plt.ylabel('height (km)')
    plt.title(r'average streamfunction over day {:2.1f} to {:2.1f}, domain size = ({:3.0f} km)$^2$'.format(t2D[0], t2D[-1], domsize))
    cb = plt.colorbar()
    cb.set_label('{:s} ({:3s})'.format(varname, vari.units.strip()))
    plt.savefig(fout + '{:s}_CRHsort_streamfnday{:2.1f}to{:2.1f}.jpg'.format(varname, t2D[0], t2D[-1]), format='jpg')
    plt.close()






      
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    