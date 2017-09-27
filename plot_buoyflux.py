from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import gc

#foutSTAT = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_STAT/'
#foutSTAT = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggr130days_64vert_ubarzero_STAT/'
foutSTAT = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggr150days_64vert_ubarzero_STAT/'
foutSTATvary = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STAT/'
#foutSTATvary = '/Users/cpatrizio/Google Drive/figures/'

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'

fout = '/Users/cpatrizio/data/SST302/'

matplotlib.rcParams.update({'font.size': 24})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'legend.fontsize': 22})
matplotlib.rcParams.update({'lines.linewidth': 3})


nc_in3D = glob.glob(fpath3D + '*512x512*day140to155*302K.nc')[0]

nc_data3D = Dataset(nc_in3D)
varis3D = nc_data3D.variables

#t3D = varis3D['time'][:]

z = varis3D['z'][:]

nz = z.size

domsizes = [768, 1536, 3072]

#start and end time in days (corresponding to first file, and last file for given domain size)
tinits = [0, 90, 110]
tends = [250, 195, 150]
colors = ['k', 'r', 'g'] 

tdouble = 90

#buoyflux = []

buoyfluxbar = []

for i, d in enumerate(domsizes):
    print 'domain size', d
    buoyfluxbar = []
    #if d == 768:
   #     fnames = glob.glob(fout + '{:d}km_buoyflux* 0to130QP.npy'.format(d))
   # elif d == 3072:
   #     fnames = glob.glob(fout + '{:d}km_buoyflux*noQP.npy'.format(d))
   # else:
    fnames = glob.glob(fout + '{:d}km_buoyflux*QP.npy'.format(d))
    for f in fnames:
        print 'loading', f
        buoyflux = np.load(f)
        ntimes = buoyflux.shape[0]
        buoyflux = buoyflux.astype('float32')
        delz = np.diff(z)
        delz2D = np.zeros((ntimes, nz-1))
        delz2D[:,:] = delz
        buoyfluxbar_zt = np.mean(np.mean(buoyflux, axis=2), axis=2)
        buoyfluxbar_int = np.sum(np.multiply(delz, buoyfluxbar_zt), axis=1)
        buoyfluxbar = np.append(buoyfluxbar, buoyfluxbar_int)
        gc.collect()
        #print 'buoyfluxbar'
        #print buoyfluxbar
    times = np.linspace(tinits[i], tends[i], len(buoyfluxbar))
    plt.figure(1)
    mask = np.isfinite(buoyfluxbar)
    plt.plot(times[mask], buoyfluxbar[mask], 'x-', color=colors[i], label='{:d} km'.format(d))
    plt.xlabel('time (days)')
    plt.ylabel('buoyancy flux (W/m$^2$)')
    plt.title('domain-mean buoyancy flux')
    gc.collect()
    #buoyfluxbar = np.append(buoyfluxbar, buoyfluxbar_int)

plt.axvline(tdouble, color='b', alpha=0.5)  
plt.legend(loc='best')
plt.savefig(foutSTATvary + 'buoyfluxNEWday0to250.png')
plt.close()


#domsize = 3072
#
##time span (days) of saved buoyancy fluxes
#tinit = 110
#tend = 140
#
#
##print 'loading file'
##np.load(fout + '{:d}km_dwdt*'.format(domsize)
#
#fnames = glob.glob(fout + '{:d}km_buoyflux*new*'.format(domsize))
#
#
#print 'domain size', domsize
#print 'stacking buoyflux arrays'
#for f in fnames:
#    print 'loading', f
#    buoyflux = np.load(f)
#    ntimes = buoyflux.shape[0]
#    delz = np.diff(z)
#    delz2D = np.zeros((ntimes, nz-1))
#    delz2D[:,:] = delz
#    buoyfluxbar_zt = np.mean(np.mean(buoyflux, axis=2), axis=2)
#    buoyfluxbar_int = np.sum(np.multiply(delz, buoyfluxbar_zt), axis=1)
#    buoyfluxbar = np.append(buoyfluxbar, buoyfluxbar_int)
#
#times = np.linspace(tinit, tend, len(buoyfluxbar))
#
#plt.figure(1)
#plt.plot(times, buoyfluxbar, 'gx-', label='3072 km')
#plt.xlabel('time (days)')
#plt.ylabel('buoyancy flux (W/m$^2$)')
#plt.title('evolution of daily averaged domain-mean buoyancy flux')
##plt.title('evolution of daily averaged domain-mean buoyancy flux, domain length = {:d} km'.format(domsize))
##plt.savefig(foutSTAT + 'buoyflux_day{:3.0f}to{:3.0f}.pdf'.format(times[0], times[-1]))
##plt.savefig(foutSTATvary + 'buoyflux.pdf')
##plt.close()
#plt.legend(loc='best')
#plt.savefig('buoyfluxnew2.pdf')
#plt.close()


    

    
