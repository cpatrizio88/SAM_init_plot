from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from SAM_init_plot.block_fns import blockave2D, blockave3D
from thermolib.constants import constants
import matplotlib.colors as colors

matplotlib.rcParams.update({'font.size': 26})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 2})
matplotlib.rcParams.update({'legend.fontsize': 22})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to170_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday90to160_64vert_ubarzero_FREQDIST/'

fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_FREQDIST/'

#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

nc_in2D1 = glob.glob(fpath2D + '*256x256*3000m*day230to250*302K.nc')[0]
nc_in3D1 = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]

nc_in2D2 = glob.glob(fpath2D +  '*512x512*3000m*day180to195*302K.nc')[0]
nc_in3D2 = glob.glob(fpath3D +  '*512x512*3000m*day180to195*302K.nc')[0]

nc_in2D3 = glob.glob(fpath2D + '*1024x1024*3000m*day170to180*302K.nc')[0]
nc_in3D3 = glob.glob(fpath3D +  '*1024x1024*3000m*day170to180*302K.nc')[0]

nc_fs = [nc_in2D1, nc_in2D2, nc_in2D3]
nc_3Dfs = [nc_in3D1, nc_in3D2, nc_in3D3]

#domsize=768
#domsize=1536
#domsize=3072

domsizes = [768, 1536, 3072]
colors=['k', 'r', 'g']

#averaging time period for 2D & 3D fields
ntave2D=24
ntave3D=4

print 'loading netCDF files'

fields = np.array([])

db=16
W_c = 0.01

varname='PW'

t=-1
nave=5
aveperiod = 24*nave
aveperiod3D = 4*nave

for i, nc_in2D in enumerate(nc_fs):
#nc_data3D = Dataset(nc_in3D) 
    nc_in3D = nc_3Dfs[i]
    print 'loading', nc_in2D
    print 'loading', nc_in3D
    domsize = domsizes[i]
    if domsize == 768:
        nave=8
    elif domsize == 1536:
        nave=8
    else:
        nave=8
    aveperiod = 24*nave
    aveperiod3D = 4*nave
    nc_data2D = Dataset(nc_in2D)
    nc_data3D = Dataset(nc_in3D)
    #varis3D = nc_data3D.variables
    varis2D = nc_data2D.variables
    varis3D = nc_data3D.variables
    x = varis2D['x'][:]
    y = varis2D['y'][:]
    #z = varis3D['z'][:]
    #p = varis3D['p'][:]
    #p = p*1e2
    W = varis3D['W'][t-aveperiod3D:t,:,:,:]
    
    Wvertbar = np.mean(np.mean(W, axis=0),axis=0)
    Wblock = blockave2D(Wvertbar, db)
    #W_c=0.5
    Wconvblock = Wblock[Wblock > W_c]
    Wdryblock = Wblock[Wblock < W_c]
    
    #W500 = varis2D['W500'][t-aveperiod:t,:,:]
    #W500_tave = np.mean(W500, axis=0)
    #W500_block = blockave2D(W500_tave, db)
    
    #nt3D = t3D.size
    #nz = z.size
    nx = x.size
    ny = y.size
    t3D = varis3D['time'][:]
    nt3D = t3D.size
    vari = varis2D[varname]
    #field = np.mean(vari[t-aveperiod:t,:,:], axis=0)
    #field = blockave2D(field, db)
    nxprime = nx / db
    nyprime = ny / db
    
    #field = field.flatten()
    #W500_flat = W500_block.flatten()
    Wconv_flat = Wconvblock.flatten()
    Wdry_flat = Wdryblock.flatten()
    nbins=100
    convbins = np.linspace(0,10,100)
    #counts, bin_edges = np.histogram(W500_flat, bins=50)
    #PWsum, bin_edges = np.histogram(W500_flat, bins=50, weights=field)
    #Wconv_counts, bin_edges = np.histogram(Wconv_flat, bins=nbins)
    Wconvsum, Wconv_bin_edges = np.histogram(Wconv_flat, bins=convbins, weights=Wconv_flat)
    Wconv_bin_c = (Wconv_bin_edges[1:] + Wconv_bin_edges[:-1])/2
    
    Wdrysum, Wdry_bin_edges = np.histogram(Wdry_flat, bins=nbins, weights=Wdry_flat)
    Wdry_bin_c = (Wdry_bin_edges[1:] + Wdry_bin_edges[:-1])/2
    
    Wconv_bars = Wconvsum/len(Wconv_flat)
    Wdry_bars = Wdrysum/len(Wdry_flat)
    
    fig = plt.figure(1)
    ax1 = fig.add_subplot(2,1,1)    
    ax1.plot(Wconv_bin_c, Wconv_bars, 'x-', color=colors[i],  label='{:d} km, day {:3.0f} to {:3.0f}'.format(domsizes[i], t3D[t-aveperiod3D], t3D[t]))
    ax1.axvline(0, color='k', alpha=0.3, linewidth=0.5)
    #axarr[0,].set_xlabel('vertically averaged W (m/s)')
    ax1.set_title(r'convective region, $w_c$ = {:2.2f} (m/s)'.format(W_c))
    ax1.set_ylabel(r'contribution to $\overline{w}_{conv}$ (m/s)')
    ax1.set_xlim((0,0.5))
    #ax.legend(loc='best')
    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(Wdry_bin_c, Wdry_bars, 'x-', color=colors[i],  label='{:d} km, day {:3.0f} to {:3.0f}'.format(domsizes[i], t3D[t-aveperiod3D], t3D[t]))
    ax2.axvline(0, color='k', alpha=0.3, linewidth=0.5)
    ax2.set_xlabel('vertically averaged W (m/s)')
    ax2.set_ylabel(r'contribution to $\overline{w}_{dry}$ (m/s)')
    #ax2.legend(loc='best')
    #plt.plot(bin_c, freqs, 'x-', color=colors[i], label='{:d} km, day {:3.0f} to {:3.0f}'.format(domsizes[i], t2D[t-aveperiod], t2D[t]))

wconv_bar = np.mean(Wconv_flat)
wdry_bar = np.mean(Wdry_flat)
    
if db == 1:
    plt.suptitle(r'$\overline{{w}}_{{conv}}$ = {:3.2f} m/s, $\overline{{w}}_{{dry}}$ = {:3.3f} m/s '.format(wconv_bar, wdry_bar))
else:
    plt.suptitle(r'block-averaging over ({:2.0f} km)$^2$, $\overline{{w}}_{{conv}}$ = {:3.2f} m/s, $\overline{{w}}_{{dry}}$ = {:3.3f} m/s '.format(db*(np.diff(x)[0])/1e3), wconv_bar, wdry_bar)
    
#if (varname == 'ZC') or (varname == 'ZE'):
#   plt.ylim(0,0.03)
#if varname == 'W500':
#   plt.xlim(0, 1)
#   plt.ylim(0, 0.1)

#plt.xlabel('vertically averaged W (m/s)')
##plt.xlim((0,0.1))
#plt.ylabel('{:s} ({:s})'.format(varname, vari.units.strip()))
#plt.ylabel('frequency')
ax1.legend(loc='best')
ax2.legend(loc='best')
plt.savefig(fout + 'convWdist_db{:d}.pdf'.format(db))
plt.close()




    
    