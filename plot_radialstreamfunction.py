from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
import SAM_init_plot.plot_misc
from SAM_init_plot.plot_misc import radprof, radprof3D
from thermolib.wsat import wsat
from thermolib.constants import constants
import SAM_init_plot.block_fns
from SAM_init_plot.block_fns import blocksort3D
from SAM_init_plot.block_fns import blockxysort2D

c = constants()

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_RADIALXSECTIONNEW/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to195_64vert_ubarzero_RADIALXSECTIONNEW/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to190_64vert_ubarzero_RADIALXSECTIONNEW/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_RADIAL/'

#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*day230to250*302K.nc')[0]

#nc_in2D = glob.glob(fpath2D + '*512x512*3000m*day180to195*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day180to195*302K.nc')[0]

nc_in2D = glob.glob(fpath2D + '*1024x1024*3000m*day170to180*302K.nc')[0]
nc_in3D = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]


#domsize=768
#domsize=1536
domsize=3072
#fraction of a day to average over.
fac=1
ndays=1

ntave2D=int(24*fac)
ntave3D=int(4*fac)

#EDIT: Time in days to look at streamfunction (t2 and t3 should correspond to the same time)
ts = np.arange(-5/fac, 0)
#ts = np.arange(30/fac, 90/fac)
ts = ts.astype(int)

#ts = [-38/fac]

#ts = [-1]



for t in ts:

    t2 = t*ntave2D
    #t2 = -1
    t3 = t*ntave3D
    #t3 = -1
    
    aveperiod2D = ntave2D
    aveperiod3D = ntave3D
    
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
    
    xx, yy = np.meshgrid(x, y)
    times = np.arange(t3D[0], np.max(t3D))
    
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
    CRH_tave = np.mean(CRH, axis=0)
    T_tave = np.mean(T, axis=0)
    w_tave = np.mean(w, axis=0)
    rhoz_tave = np.mean(rhoz, axis=0)
    
    
    #number of grid points in a block
    db=16

    PW_tave = np.mean(PW, axis=0)
    PWxy_sorted = blockxysort2D(PW_tave, xx, yy, db)
    PWsort = PWxy_sorted.keys()
    PWxycoords = PWxy_sorted.values()

    mcenter = PWxycoords[-1]
    
    binwidth=5e3
    bins = np.arange(0, domsize*1e3, binwidth)
    
    
    #wrsort = np.zeros(z.size, bins.size)
    psi = np.zeros((z.size, bins.size-1))
    chi = np.zeros((z.size, bins.size-1))
    
    #calculate the stream function in cylindrical coords 
    for i, zlev in enumerate(z):
        wrbins, wmeans = radprof(w_tave[i,:,:], xx, yy, mcenter, bins)
        delr = np.diff(wrbins)
        wrbin_centers = (wrbins[:-1] + wrbins[1:])/2.
        chi[i,:] = rhoz_tave[i]*np.cumsum(np.multiply(wrbin_centers*delr, wmeans))
        psi[i,:] = chi[i,:]/wrbin_centers
        
        #for j, r in enumerate(wrbin_centers[:-1]):
        #  r1 = wrbin_centers[j+1]
        #  r0 = wrbin_centers[j]
          #psi[i,j+1] = (r0/r1)*(psi[i, j] + rhoz_tave[i]*np.cumsum(np.multiply(delr, wmeans))
        #  psi[i,j+1] = (r0/r1)*(psi[i, j] + rhoz_tave[i]*wmeans[j]*delr[j])

        
    
    #EDIT: variable to overlay with streamfunction
    varname = 'QV'
    vari = varis3D[varname]
    field = vari[t3-aveperiod3D:t3,:,:,:]
    field_tave = np.mean(field, axis=0)
    
    zedges = np.zeros(z.shape)
    zedges[1:] = (z[1:] + z[:-1])/2.
    
    bins3D = [bins, zedges]
    
    fieldrbins, zedges, fieldmeans = radprof3D(field_tave, xx, yy, z, mcenter, bins3D)
    fieldrbin_centers = (fieldrbins[:-1] + fieldrbins[1:])/2.
    zcenters = (zedges[:-1] + zedges[1:])/2.
    
    rrpsi, zzpsi = np.meshgrid(wrbin_centers, z)
    rr, zz = np.meshgrid(fieldrbin_centers, zcenters)
    
    if varname=='QV':
       #cvals = np.arange(0, 19)
       vmin=0
       vmax=19
    elif varname=='QN':
       vmin=0
       vmax=0.2
    else:
       cvals=20
    
    fieldmeans = np.transpose(fieldmeans)
    plt.figure()
    ax=plt.gcf().gca()
    plt.contour(rrpsi/(domsize*1e3), zzpsi/1e3, psi, 15, linewidths=1.5, colors='k', zorder=1)
    plt.pcolormesh(rr/(domsize*1e3), zz/1e3, fieldmeans, vmin=vmin, vmax=vmax, cmap=cm.RdYlBu_r, zorder=0)
    ax.set_ylim((0,16))
    plt.xlabel(r'$\hat{r}$')
    #plt.xlabel('radial distance (km)')
    plt.ylabel('z (km)')
    plt.title(r'streamfunction, day {:2.0f} to {:2.0f} average, domain size = ({:3.0f} km)$^2$'.format(t2D[0], t2D[-1], domsize))
    cb = plt.colorbar()
    ax.set_xlim([0, 1./np.sqrt(2)])
    #ax.set_xlim([0, (1./np.sqrt(2))*domsize])
    cb.set_label('{:s} ({:3s})'.format(varname, vari.units.strip()))
    plt.savefig(fout + '{:s}_radial_streamfnday{:2.1f}to{:2.1f}.jpg'.format(varname, t2D[0], t2D[-1]), format='jpg')
    plt.close()

    
    
    
        
        
        
    
        
        
        
    
    
    
   