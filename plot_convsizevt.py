from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from SAM_init_plot.block_fns import blockave2D, blockave3D
from thermolib.constants import constants
import matplotlib.colors as colors

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to170_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday90to160_64vert_ubarzero_FREQDIST/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STAT/'

#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

#nc_in2D = glob.glob(fpath2D + '*1024x1024*3000m*day090to130*302K.nc')[0]

nc_fs1 = glob.glob(fpath2D + '*256x256*3000m*302K.nc')
nc_fs2 = glob.glob(fpath2D + '*512x512*3000m*302K.nc')
nc_fs3 = glob.glob(fpath2D + '*1024x1024*3000m*302K.nc')
domsize=768
#domsize=1536
#domsize=3072

nc_fs = [nc_fs1, nc_fs2, nc_fs3]
domsizes=[768, 1536, 3072]
colors = ['k', 'r', 'g']


#z = varis3D['z'][:]
#p = varis3D['p'][:]
#p = p*1e2

tstarts = [0,90,90]

for i, nc_f in enumerate(nc_fs):
        
    print 'domsize', domsizes[i]
    
    print 'loading netCDF files'
    
    convfracs = np.array([])
    
    db=1
    W_c = 0.5
    
    times = np.array([])
    
    varname='W500'

    for nc_in2D in nc_f:
    #nc_data3D = Dataset(nc_in3D) 
        
        print 'loading', nc_in2D
        
        nc_data2D = Dataset(nc_in2D)
        varis2D = nc_data2D.variables
        
        x = varis2D['x'][:]
        y = varis2D['y'][:]
        
        #nt3D = t3D.size
        #nz = z.size
        nx = x.size
        ny = y.size
        
        nxprime = nx / db
        nyprime = ny / db
        
        #averaging time period for 2D & 3D fields
        ntave2D=24
        ntave3D=4
        nc_data2D = Dataset(nc_in2D)
        #varis3D = nc_data3D.variables
        varis2D = nc_data2D.variables
        t2D = varis2D['time'][:]
        nt2D = t2D.size
        field = varis2D[varname][:]
        
        #ntrunc = nt2D%ntave2D
        #ntrunc=0
        #field = field[ntrunc:,:,:]
        units = varis2D[varname].units
        #field_tmp = field.reshape(nt2D/ntave2D, ntave2D, nx, ny)
        #field_tave = np.mean(field_tmp, axis=1)
        if db == 1:
            field_blockave = field
        else:
            field_blockave = blockave3D(field, db)
        nt = field_blockave.shape[0]
        convfrac = np.zeros(nt)
        #subsfrac = np.zeros(nt)
        for k, t in enumerate(t2D):
            field_t = field_blockave[k,:,:]
            convfrac[k] = len(field_t[field_t[:,:] > W_c])/(1.*nxprime*nyprime)
            #subsfrac[i] = 1 - convfrac[i]
        
        #        buoyflux = np.vstack((buoyflux, buoyflux_temp)) if buoyflux.size else buoyflux_temp
        convfracs = np.concatenate((convfracs, convfrac)) if convfracs.size else convfrac
        times = np.concatenate((times, t2D)) if times.size else t2D
        
    #ntimes = convfracs.shape[0]
    
    print 'plotting'
    
    #calculate 2D (CRH rank, time) map of CRH frequency  ###
    #times = np.arange(tstarts[i], tstarts[i]+ntimes)

    w = 24*5
    
    convfracs_smooth = moving_average(convfracs, w)
    
    plt.figure(1)
    plt.plot(times, convfracs, color=colors[i], alpha=0.5, label='{:d} km'.format(domsizes[i]))
    plt.plot(times[w/2:-w/2+1], convfracs_smooth, color=colors[i], linewidth=1.5)
    #plt.plot(times, 1-convfracs, 'k--')
    plt.xlabel('time (days)')
    plt.ylabel('convective fractional area')
    plt.ylim((0, 0.10))
    if db == 1:
        plt.title('convective fractional area (W500 > W$_c$ = {:3.2f} m/s)'.format(W_c))
    else:
        plt.title('convective fractional area (W500 > W$_c$ = {:3.2f} m/s), average over ({:2.0f} km)$^2$ blocks '.format(W_c, db*np.diff(x)[0]/1e3))
    plt.savefig(fout + 'fracconv250days_db{:d}.pdf'.format(db))
  

plt.legend()
plt.savefig(fout + 'fracconv250days_db{:d}.pdf'.format(db))
plt.close()
        

    