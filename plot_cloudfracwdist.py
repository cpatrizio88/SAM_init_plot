from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker
import matplotlib

from SAM_init_plot.block_fns import blocksum3D, blockave3D, blockave2D
from thermolib.constants import constants


matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 2})

plt.style.use('seaborn-white')
fpathSTAT = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
fpath2D =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fpath3D = '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
foutdata = '/Users/cpatrizio/data/SST302/'
#fout =  '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_MAPSNEW/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to140_64vert_ubarzero_MAPSNEW/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to150_64vert_ubarzero_MAPSNEW/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_CLOUDFRAC/'

nc_STAT = glob.glob(fpathSTAT + '*256x256*3000m*250days*302K.nc')[0]
nc_in2D1 = glob.glob(fpath2D + '*256x256*3000m*day230to250*302K.nc')[0]
nc_in3D1 = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]
 
#nc_STAT = glob.glob(fpathSTAT + '*512x512*3000m*180days*302K.nc')[0]
nc_in2D2 = glob.glob(fpath2D + '*512x512*3000m*day170to180*302K.nc')[0]
nc_in3D2 = glob.glob(fpath3D + '*512x512*3000m*day170to180*302K.nc')[0]

#nc_STAT = glob.glob(fpathSTAT + '*1024x1024*3000m*170days*302K.nc')[0]
nc_in2D3 = glob.glob(fpath2D + '*1024x1024*3000m*day170to180*302K.nc')[0]
nc_in3D3 = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]

domsize=768
#domsize=1536
#domsize=3072

domsizes = [768, 1536, 3072]
colors = ['k', 'r', 'g']
nc_2Ds = [nc_in2D1, nc_in2D2, nc_in2D3]
nc_3Ds = [nc_in3D1, nc_in3D2, nc_in3D3]

#varname = 'W500'
#vari = varis[varname]
#field = vari[:]

for i, domsize in enumerate(domsizes):
    
    print 'domsize', domsize

    nc_data2D= Dataset(nc_2Ds[i])
    nc_dataSTAT = Dataset(nc_STAT)
    nc_data3D = Dataset(nc_3Ds[i])
    varis2D = nc_data2D.variables
    varis3D = nc_data3D.variables
    varisSTAT = nc_dataSTAT.variables
    times2D = varis2D['time'][:]
    
    #times3D = varis3D['time'][:]
    
    x = varis2D['x'][:]
    y = varis2D['y'][:]
    
    p = varis3D['p'][:]
    z = varis3D['z'][:]
    
    nx = x.size
    ny = y.size
    
    #fac=8
    
    nave=9
    
    nave2D = 24
    nave3D = 4
    
    aveperiod2D=nave2D*nave
    aveperiod3D=nave3D*nave
    
    t = -1
    
    db=16
    
    nxprime = nx // db
    
    varname = 'W500'
    
    W = varis2D[varname][:]
    
    #W = varis3D['W'][t-nave3D:t,:,:,:]
    
    W = varis2D[varname][t-aveperiod2D:t,:,:]
    W_tave = np.mean(W, axis=0)
    
    W_blockave = blockave2D(W_tave, db)
    
    
    #put code to select W at specific level
    #W_blockave = blockave2D(W_tavez, db)
    
    
    QN = varis3D['QN'][:]
    
    QN_tave = np.mean(QN[t-aveperiod3D:t,:,:,:], axis=0)
    
    cloudfreq = np.sign(QN_tave)
    
    cloudcount = blocksum3D(cloudfreq, db)
    
    cloudfrac = cloudcount/(db*db)
    
    plow_ti = np.where(p > 680)[0][-1]
    plow_bi = np.where(p < 1000)[0][0]
    
    phigh_bi = np.where(p < 300)[0][0]
    phigh_ti = np.where(p < 100)[0][0]
    
    #lowcloudfrac = np.mean(cloudfrac[plow_bi:plow_ti,:,:], axis=0)
    #midcloudfrac = np.mean(cloudfrac[plow_ti:phigh_bi,:,:], axis=0)
    #highcloudfrac = np.mean(cloudfrac[phigh_bi:phigh_ti,:,:], axis=0)
    
        
    lowcloudfrac = np.mean(cloudfrac[plow_bi:plow_ti,:,:], axis=0)
    midcloudfrac = np.mean(cloudfrac[plow_ti:phigh_bi,:,:], axis=0)
    highcloudfrac = np.mean(cloudfrac[phigh_bi:phigh_ti,:,:], axis=0)
    
    wflat = W_blockave.flatten()
    lowcloudfracflat = 100*lowcloudfrac.flatten()
    midcloudfracflat = 100*midcloudfrac.flatten()
    highcloudfracflat = 100*highcloudfrac.flatten()
    
    nxprime = nx // db
    
    wbins = np.linspace(0, 1, 100)
    
    print 'binning vertical velocity'
    
    nwhist, wedges = np.histogram(wflat, bins=wbins)
    whist, wedges = np.histogram(wflat, bins=wbins, weights=np.zeros_like(wflat) + 1. / (nxprime*nxprime))
    
    print 'binning cloud fraction'
    lowcldhist, wedges = np.histogram(wflat, bins=wbins, weights=lowcloudfracflat)
    midcldhist, wedges = np.histogram(wflat, bins=wbins, weights=midcloudfracflat)
    highcldhist, wedges = np.histogram(wflat, bins=wbins, weights=highcloudfracflat)
    
    lowcldhist = lowcldhist/nwhist
    midcldhist = midcldhist/nwhist
    highcldhist = highcldhist/nwhist
    
    wcenters = (wedges[:-1] + wedges[1:])/2.
    
    indices = np.isfinite(lowcldhist)
    
    cloudfrac_z = np.mean(np.mean(cloudfrac, axis=1), axis=1)
    
    plt.figure(1)
    plt.plot(cloudfrac_z, z, color = colors[i])
    plt.axhline(z[plow_bi], alpha=0.5)
    plt.axhline(z[plow_ti], alpha=0.5)
    plt.axhline(z[phigh_bi], alpha=0.5)
    plt.axhline(z[phigh_ti], alpha=0.5)
    plt.show()
    
    
    fig=plt.figure(2)
    plt.suptitle('W500 frequency distribution (top panel) and cloud fraction (bottom panels)')
    ax1 = fig.add_subplot(4,1,1)
    ax1.plot(wcenters[indices], whist[indices], '-', color=colors[i], linewidth=3, label='{:d} km, day {:3.0f} to day {:3.0f} average'.format(domsize, times2D[t-aveperiod2D], times2D[t]))
    ax1.set_ylabel('frequency')
    plt.axvline(0, color='k', alpha=0.5)
    ax2 = fig.add_subplot(4,1,2)
    ax2.plot(wcenters[indices], lowcldhist[indices], 'x-', color=colors[i], label='{:d} km, day {:3.0f} to day {:3.0f} average'.format(domsize, times2D[t-aveperiod2D], times2D[t]))
    ax2.set_ylabel('low cloud fraction (%)')
    plt.axvline(0, color='k', alpha=0.5)
    ax3 = fig.add_subplot(4,1,3)
    ax3.plot(wcenters[indices], midcldhist[indices], 'o-', color=colors[i], label='{:d} km, day {:3.0f} to day {:3.0f} average'.format(domsize, times2D[t-aveperiod2D], times2D[t]))
    ax3.set_ylabel('mid cloud fraction (%)')
    plt.axvline(0, color='k', alpha=0.5)
    ax4 = fig.add_subplot(4,1,4)
    ax4.plot(wcenters[indices], highcldhist[indices], '^-', color=colors[i], label='{:d} km, day {:3.0f} to day {:3.0f} average'.format(domsize, times2D[t-aveperiod2D], times2D[t]))
    ax4.set_ylabel('high cloud fraction (%)')
    plt.axvline(0, color='k', alpha=0.5)
    ax2.set_ylim(0,100)
    ax3.set_ylim(0,100)
    ax4.set_ylim(0,100)
    ax4.set_xlabel('W500 (m/s)')
    

plt.legend(loc='best')
plt.savefig(fout + 'W500vscloudfrac.pdf')



















#block average W500, find the vertically averaged cloud fraction at different levels for each block ?






