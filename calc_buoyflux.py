# -*- coding: utf-8 -*-
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
import SAM_init_plot.block_fns
import SAM_init_plot.misc_fns
from thermolib.constants import constants
from SAM_init_plot.block_fns import blockave2D, blockave3D, blockxysort2D
from SAM_init_plot.misc_fns import radprof3D
import gc

c=constants()

eps = c.Rd/c.Rv
epsv = (1-eps)/eps

matplotlib.rcParams.update({'font.size': 28})
matplotlib.rcParams.update({'figure.figsize': (18, 10)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 22})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/data/SST302/'
#foutSTAT = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_STAT/'
#foutSTAT = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggr130days_64vert_ubarzero_STAT/'
foutSTAT = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggr150days_64vert_ubarzero_STAT/'

nc_in3D = glob.glob(fpath3D + '*256x256*3000m*day130days*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day170to180*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*1024x1024*3000m*day140to150*302K.nc')[0]

domsize=768
#domsize=1536
#domsize=3072

nc_data3D = Dataset(nc_in3D)
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




dz=0.5*(z[0] + z[1])
nzm=nz-1

#coords of midpoint
adz = np.zeros(nzm)
for k in np.arange(1, nzm-1):
   adz[k] = 0.5*(z[k+1]-z[k-1])/dz
adz[0] = 1
adz[-1] = (z[-1] - z[-2])/dz

#betu = adz[:-1]/(adz[1:] + adz[:-1])
#betd = adz[1:]/(adz[1:] + adz[:-1])

#ts = np.arange(1, 10)
#ntimes = len(ts)

#tstart = 0*ntave3D
#tend = 100*ntave3D
#tend= -1
tstart = 0
tend= -1

t3D = t3D[tstart:tend]
#t3D = t3D[:11*ntave3D]
#t3D = t3D[11*ntave3D:16*ntave3D]
#t3D = t3D[11*ntave3D:]
#t3D = t3D[16*ntave3D:]
ntimes = len(t3D)

#LOOP OVER TIME HERE?
dwdt = np.zeros((ntimes, nz, nx, ny))
buoyflux = np.zeros((ntimes, nz, nx, ny))
#Tprimeu = np.zeros(buoyflux.shape)
#Tprimed = np.zeros(buoyflux.shape)

print 'domsize', domsize


for i, t in enumerate(t3D):
    #t3=t*ntave3D
    t3=tstart + i
    print 'calculating dw/dt at day {:3.2f}'.format(t)
    #print 'calculating dw/dt between day {:3.2f} and {:3.2f}'.format(t3D[t3-ntave3D], t3D[t3])
    qv = varis3D['QV'][t3,:,:,:]
    qv = qv.astype('float32')
    qv = qv*1e-3
    qn = varis3D['QN'][t3,:,:,:]
    qn = qn.astype('float32')
    qn = qn*1e-3
    #IGNORE PRECIP LOADING BECAUSE LARGE DOMAIN IS MISSING QP FIELD.
    #if domsize == 3072:
    #    #temporary fix for missing QP -> ignore precip loading. 
    #    qp = np.zeros(qn.shape)
    #else:
    qp = varis3D['QP'][t3,:,:,:]
    qp = qp*1e-3
    
    qp = qp.astype('float32')    
    T = varis3D['TABS'][t3,:,:,:]
    T = T.astype('float32')
    w = varis3D['W'][t3,:,:,:]
    w = w.astype('float32')
    
    rhoz = p/(c.Rd*np.mean(np.mean(T, axis=2), axis=1))
    T_tave = T
    qv_tave = qv
    qn_tave = qn
    qp_tave = qp
    w_tave = w
    rhoz_tave = rhoz
    #qv_tave = np.mean(qv, axis=0)
    #qn_tave = np.mean(qn, axis=0)
    #qp_tave = np.mean(qn, axis=0)
    #T_tave = np.mean(T, axis=0)
    #w_tave = np.mean(w, axis=0)
    #rhoz_tave = np.mean(rhoz, axis=0)
    T0_tave = np.mean(np.mean(T_tave, axis=2), axis=1)
    qv0_tave = np.mean(np.mean(qv_tave, axis=2), axis=1)
    qn0_tave = np.mean(np.mean(qn_tave, axis=2), axis=1)
    qp0_tave = np.mean(np.mean(qp_tave, axis=2), axis=1)
    #theta0_tave = T0_tave*(c.p0/p)**(c.Rd/c.cp)
    #thetav0_tave = theta0_tave*(1+epsv*qv0_tave)
    bet = c.g/T0_tave
    
    for k in np.arange(1, nzm):
        #print 'level', k
        kb = k-1
        betu = adz[kb]/(adz[k] + adz[kb])
        betd = adz[k]/(adz[k] + adz[kb])
        qvprimeu = qv_tave[k,:,:] - qv0_tave[k]
        qnprimeu = qn_tave[k,:,:] - qn0_tave[k]
        qpprimeu = qp_tave[k,:,:] - qp0_tave[k]
        Tprimeu = T_tave[k,:,:] - T0_tave[k]
        qvprimed = qv_tave[kb,:,:] - qv0_tave[kb]
        qnprimed = qn_tave[kb,:,:] - qn0_tave[kb]
        qpprimed = qp_tave[kb,:,:] - qp0_tave[kb]
        Tprimed = T_tave[kb,:,:] - T0_tave[kb]
        
        dwdt[i, k,:,:] = bet[k]*betu* \
                    (T0_tave[k]*(epsv*qvprimeu -(qnprimeu + qpprimeu)) + \
                    Tprimeu*(1 + epsv*qv0_tave[k] - qn0_tave[k] - qp0_tave[k])) + \
                    bet[kb]*betd* \
                    (T0_tave[kb]*(epsv*qvprimed -(qnprimed + qpprimed)) + \
                    Tprimed*(1 + epsv*qv0_tave[kb] - qn0_tave[kb] - qp0_tave[kb]))
        
        buoyflux[i, k,:,:] = rhoz_tave[k]*w_tave[k,:,:]*dwdt[i, k,:,:]
        dwdt[i,k,:,:] = rhoz_tave[k]*dwdt[i,k,:,:]
        gc.collect()
    gc.collect()
 
#eliminate first row of zeros (cant compute d/dz at surface)                               
dwdt = dwdt[:,1:,:,:]
buoyflux = buoyflux[:,1:,:,:]

trunc=ntimes%ntave3D

dwdt = dwdt[trunc:,:,:,:]
#buoyflux=buoyflux[trunc:,:,:,:]

times = np.arange(int(t3D[0]), int(t3D[-1]))

#dwdtnew = np.zeros((times.size, nz-1, nx, ny))
buoyfluxnew = np.zeros((times.size, nz-1, nx, ny))

print 'calculating daily averages'
for ti, t in enumerate(times):
    print 'day', t
    #dwdtnew[ti,:,:,:] = np.mean(dwdt[ti*ntave3D:(ti+1)*ntave3D,:,:,:], axis=0)
    buoyfluxnew[ti,:,:,:] = np.mean(buoyflux[ti*ntave3D:(ti+1)*ntave3D,:,:,:], axis=0)
    
#. It's as simple as f.create_dataset('name', data=x) where x is your numpy array and f is the open hdf file. 
# Doing the same thing in pytables is possible, but considerably more difficult. –
    

#dwdt = dwdt.reshape(ntimes/ntave3D, ntave3D, nz-1, nx, ny)
#dwdt = np.mean(dwdt, axis=1)

#buoyflux = buoyflux.reshape(ntimes/ntave3D, ntave3D, nz-1, nx, ny)
#buoyflux = np.mean(buoyflux, axis=1)

print 'saving buoyancy file'

#np.save(fout + '{:d}km_buoy_day{:3.0f}to{:3.0f}'.format(domsize, t3D[0], t3D[-1]), dwdtnew)

#print 'saving buoyancy flux file'
#np.save(fout + '{:d}km_dwdt_day{:3.0f}to{:3.0f}QP'.format(domsize, t3D[0], t3D[-1]), dwdtnew)
np.save(fout + '{:d}km_buoyflux_day{:3.0f}to{:3.0f}'.format(domsize, t3D[0], t3D[-1]), buoyfluxnew)

delz = np.diff(z)
delz2D = np.zeros((ntimes, nz-1))
delz2D[:,:] = delz

#buoyfluxbar_zt = np.mean(np.mean(buoyflux, axis=2), axis=2)

#buoyfluxbar = np.sum(np.multiply(delz, buoyfluxbar_zt), axis=1)

buoybar_zt = np.mean(np.mean(dwdt, axis=2), axis=2) #has units N/m^3

buoybar = np.sum(np.multiply(delz, buoybar_zt), axis=1) #has units N/m^2 = Pa

#print 'mean column-integrated buoyancy flux = {:3.2f} (W/m^2)'.format(buoyfluxbar)

times = t3D

times = times[trunc:]

#plt.figure(1)
#plt.plot(times, buoyfluxbar, 'kx-')
#plt.ylabel(r'column-integrated buoyancy flux (W/m$^2$)')
#plt.xlabel('t (days)')
#plt.title(r'evolution of domain-mean buoyancy flux, domain size ({:d} km)$^2$'.format(domsize)) 
#plt.savefig(foutSTAT + 'buoyflux_day{:3.0f}to{:3.0f}QP.pdf'.format(t3D[0], t3D[-1]))
#plt.close()

#plt.figure(2)
#plt.plot(times, buoybar, 'kx-')
#plt.ylabel(r'column-integrated buoyancy  (Pa)')
#plt.xlabel('t (days)')
#plt.title(r'evolution of domain-mean buoyancy, domain size ({:d} km)$^2$'.format(domsize)) 
#plt.savefig(foutSTAT + 'buoy_day{:3.0f}to{:3.0f}.pdf'.format(t3D[0], t3D[-1]))
#plt.close()




                  
                  
                
    
    
    







