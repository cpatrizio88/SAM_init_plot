from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
import SAM_init_plot.block_fns
import SAM_init_plot.misc_fns
from thermolib.constants import constants
import gc

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/data/SST302/'
foutSTAT = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STAT/'

domsizes=[768, 1536, 3072]
colors = ['k', 'r', 'g']
tinits = [0, 90, 110]
tends = [250, 140, 150]

for i, d in enumerate(domsizes):
    
    KE_bar = np.array([])
    times = np.array([])
    
    print 'domain size', d
    fnames = glob.glob(fout + '{:d}km_KEbar*'.format(d))
    for fname in fnames:
        KE_bar_temp = np.load(fname)
        KE_bar = np.append(KE_bar, KE_bar_temp)
    
    ntimes = KE_bar.shape[0]  
    times = np.linspace(tinits[i], tends[i], ntimes)              
    plt.figure(1)
    plt.plot(times, KE_bar, color=colors[i], label='{:d} km'.format(d))
    plt.xlabel('time (days)')
    plt.ylabel(r'mean kinetic energy (J/m$^2$)')
    plt.title(r'evolution of domain-mean kinetic energy')
    plt.savefig(foutSTAT + 'KE_day0to250.pdf')
    
        
tdouble=90

plt.axvline(tdouble, color='b', alpha=0.5)  
plt.legend(loc='best')
plt.savefig(foutSTAT + 'KE_day0to250.pdf')
plt.close()
    

    