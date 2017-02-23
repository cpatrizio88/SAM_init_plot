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
fpathSTAT =  '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to170_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday90to160_64vert_ubarzero_FREQDIST/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STAT/'

#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

#nc_in2D = glob.glob(fpath2D + '*1024x1024*3000m*day090to130*302K.nc')[0]

nc_f2D1 = glob.glob(fpath2D + '*256x256*3000m*302K.nc')
nc_f2D2 = glob.glob(fpath2D + '*512x512*3000m*302K.nc')
nc_f2D3 = glob.glob(fpath2D + '*1024x1024*3000m*302K.nc')

nc_fSTAT1 = glob.glob(fpathSTAT + '*256x256*3000m*250days_302K.nc')[0]
nc_fSTAT2 = glob.glob(fpathSTAT + '*512x512*3000m*195days_302K.nc')[0]
nc_fSTAT3 = glob.glob(fpathSTAT + '*1024x1024*3000m*200days_302K.nc')[0]

#domsize=768
#domsize=1536
#domsize=3072

nc_f2Ds = [nc_f2D1, nc_f2D2, nc_f2D3]
nc_fSTATs = [nc_fSTAT1, nc_fSTAT2, nc_fSTAT3]
domsizes=[768, 1536, 3072]
colors = ['k', 'r', 'g']

#nc_f2Ds = [nc_f2D1]
#nc_fSTATs = [nc_fSTAT1]
#domsizes=[768]

tstarts = [0,90,90]

                
db = 1
W_crit = 0.5

for i, nc_flist in enumerate(nc_f2Ds):
        
    print 'domsize', domsizes[i]
    convvals = []
    dryvals = []
    times = np.array([])
    
    for nc_f in nc_flist:
    

    
    
        varname='LWNT'
    
    
        print 'loading', nc_f
        
        #averaging time period for 2D & 3D fields
        ntave2D=24
        ntave3D=4
        nc_data2D = Dataset(nc_f)
        varis2D = nc_data2D.variables
        t2D = varis2D['time'][:]
        nt2D = t2D.size
        vari = varis2D[varname]
        field = vari[:]
        fieldblock = blockave3D(field, db)
        
        W500 = varis2D['W500'][:]
        
 
        W500block = blockave3D(W500, db)
        
        nt = W500block.shape[0]
    
   #convective region and dry region time series
        print 'calculating convective region + dry region averages'
        for k,t in enumerate(t2D):

            W500block_t = W500block[k,:,:]
            conv_points = W500block_t >= W_crit
            dry_points = W500block_t < W_crit
            
            conv_t = np.mean(fieldblock[t, conv_points])
            dry_t = np.mean(fieldblock[t, dry_points])
            convvals = convvals + [conv_t]
            dryvals = dryvals + [dry_t]
            
            

            
        
        times = np.concatenate((times, t2D)) if times.size else t2D
        
    convvals = np.array(convvals)
    dryvals = np.array(dryvals)    
    
    units = vari.units
    
    print 'plotting'
    w = 24*5
    
    x = varis2D['x'][:]
    
    convvals[np.isnan(convvals)] = np.mean(convvals[~np.isnan(convvals)])
    dryvals[np.isnan(dryvals)] = np.mean(dryvals[~np.isnan(dryvals)])
    
    convvals_smooth = moving_average(convvals, w)
    dryvals_smooth = moving_average(dryvals, w)
    
    plt.figure(1)
    plt.plot(times, convvals, color=colors[i], alpha=0.5, label='{:d} km'.format(domsizes[i]))
    plt.plot(times[w/2:-w/2+1], convvals_smooth, color=colors[i], linewidth=1.5)
    plt.xlabel('time (days)')
    plt.ylabel('{:s} ({:s})'.format(varname, units))
    if db == 1:
        plt.title('average {:s} over convective region (W$_c$ > {:3.2f} m/s)'.format(varname, W_crit))
    else:
        plt.title('average {:s} over convective region (W$_c$ > {:3.2f} m/s), block-averaging over ({:2.0f} km)$^2$'.format(varname, W_crit, db*(np.diff(x)[0])/1e3))
    plt.savefig(fout + 'conv{:s}250days_db{:d}.pdf'.format(varname, db))
    
    plt.figure(2)
    plt.plot(times, dryvals, color=colors[i], alpha=0.5, label='{:d} km'.format(domsizes[i]))
    plt.plot(times[w/2:-w/2+1], dryvals_smooth, color=colors[i], linewidth=1.5)
    plt.xlabel('time (days)')
    plt.ylabel('{:s} ({:s})'.format(varname, units))
    plt.title('average {:s} over convection-free region (W$_c$ < {:3.2f} m/s), block-averaging over ({:2.0f} km)$^2$'.format(varname, W_crit, db*(np.diff(x)[0])/1e3))
    #plt.ylim((-.05,0))
    plt.savefig(fout + 'dry{:s}250days_db{:d}.pdf'.format(varname, db))
  
plt.figure(1)
plt.legend()
plt.savefig(fout + 'conv{:s}250days_db{:d}.pdf'.format(varname, db))
plt.close()

plt.figure(2)
plt.legend()
plt.savefig(fout + 'dry{:s}250days_db{:d}.pdf'.format(varname, db))
plt.close()

        

    