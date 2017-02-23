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

c = constants()

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/data/SST302/'
foutSTAT = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STAT/'
#foutSTAT = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggr130days_64vert_ubarzero_STAT/'
#foutSTAT = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggr150days_64vert_ubarzero_STAT/'

#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day130to140*302K.nc')[0]
nc_in3D = glob.glob(fpath3D + '*1024x1024*3000m*day110to120*302K.nc')[0]

#domsize=768
#domsize=1536
domsize=3072

domsizes=[1024]
colors = ['k', 'r', 'g']

for i, d in enumerate(domsizes):
    
    KE_bar = np.array([])
    times = np.array([])
    
    print 'domain size', d*3
    
    nc_in3D = glob.glob(fpath3D + '*{:d}*302K.nc'.format(d))
    
    print 'files: ', nc_in3D
    
    for f in nc_in3D:
        
        print 'loading', f
   
        nc_data3D = Dataset(f)
        #nc_data2D = Dataset(nc_in2D)
        varis3D = nc_data3D.variables
        
        t3D = varis3D['time'][:]
        #t2D = varis2D['time'][:]
        x = varis3D['x'][:]
        y = varis3D['y'][:]
        z = varis3D['z'][:]
        p = varis3D['p'][:]
        p = p*1e2
        
        ntave3D=4
        nt3D = t3D.size
        nz = z.size
        nx = x.size
        ny = y.size
        
        #tstart = 0
        #tend= -1
        
        #t3D = t3D[tstart:tend]
        #t3D = t3D[:11*ntave3D]
        #t3D = t3D[11*ntave3D:16*ntave3D]
        #t3D = t3D[11*ntave3D:]
        #t3D = t3D[16*ntave3D:]
        ntimes = len(t3D)
        
        KE = np.zeros((ntimes, nz, nx, ny))
        
        p3D = np.zeros((nz, nx, ny))
        p3D = p3D.T
        p3D[:,:,:] = p
        p3D = p3D.T
        #fnames = glob.glob(fout + '{:d}km_buoyflux*QP.npy'.format(d))
    

        for ti, t in enumerate(t3D):
        #t3=t*ntave3D
            #t3=tstart + i
            print 'calculating KE at day {:3.2f}'.format(t)
            
            W = varis3D['W'][ti,:,:,:]
            U = varis3D['U'][ti,:,:,:]
            V = varis3D['V'][ti,:,:,:]
            
            T = varis3D['TABS'][ti,:,:,:]
            
        
            
            rho = p3D/(c.Rd*T)
        
        
            
            KE[ti,:,:,:] = 0.5*rho*(W**2 + U**2 + V**2)
            #KE[ti,:,:,:] = 0.5*rho*W**2
            gc.collect()
    
        KE_barzt = np.mean(np.mean(KE, axis=3), axis=2)
        delz = np.diff(z)
        delz2D = np.zeros((ntimes, nz-1))
        delz2D[:,:] = delz
        KE_bar_temp = np.sum(np.multiply(delz, KE_barzt[:,:-1]), axis=1)
        KE_bar = np.append(KE_bar, KE_bar_temp)
        times = np.append(times, t3D)
        gc.collect()
        
    #print 'domsize', d*3
    #print 'length 
    
    np.save(fout + '{:d}km_KEbar_day{:3.0f}to{:3.0f}'.format(d*3, times[0], times[-1]), KE_bar)
    
#    plt.figure(1)
#    plt.plot(times, KE_bar, color=colors[i], label='{:d} km'.format(3*d))
#    plt.xlabel('time (days)')
#    plt.ylabel(r'mean vertical kinetic energy (J/m$^2$)')
#    plt.title(r'evolution of domain-mean vertical kinetic energy')
#    plt.savefig(foutSTAT + 'wKE_day0to250.pdf')
#    
#tdouble=90
#
#plt.axvline(tdouble, color='b', alpha=0.5)  
#plt.legend(loc='best')
#plt.savefig(foutSTAT + 'KE_day0to250.pdf')
#plt.close()
    

    
    
    