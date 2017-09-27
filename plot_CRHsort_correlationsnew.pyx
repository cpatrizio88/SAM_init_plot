from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from thermolib.wsat import wsat
from thermolib.constants import constants
import SAM_init_plot.block_fns
from SAM_init_plot.block_fns import blockave2D, blocksort2D, blocksort3D, vertint
import collections

import pyximport; pyximport.install(pyimport = True)
cimport cython.view
from cython.view cimport array as cvarray


c = constants()

p_s = 1015
T_s = 302


def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (22, 14)})
matplotlib.rcParams.update({'lines.linewidth': 2})
matplotlib.rcParams.update({'legend.fontsize': 22})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_CRHSORT/'
#fout = '/Users/cpatrizio/figures/SST302/1536km_SAM_aggr140days_64vert_ubarzero_CRHSORT'


nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

#nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day090to130*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*512x512*3000m*day090to140*302K.nc')[0]

nc_data3D = Dataset(nc_in3D)
nc_data2D = Dataset(nc_in2D)
varis3D = nc_data3D.variables
varis2D = nc_data2D.variables

t3D = varis3D['time'][:]
t2D = varis2D['time'][:]
x = varis3D['x'][:]
y = varis3D['y'][:]
z = varis3D['z'][:]
p = varis3D['p'][:]
p = p*1e2
#2D variables
print 'loading 2D variables'
#PW = varis2D['PW'][:]
cdef float[:,:,:,:] PW = varis2D['PW'][:]
IWP = varis2D['IWP'][:]
LHF = varis2D['LHF'][:]
SHF = varis2D['SHF'][:]
SWVP = varis2D['SWVP'][:]
SWNT = varis2D['SWNT'][:]
SWNS = varis2D['SWNS'][:]
LWNT = varis2D['LWNT'][:]
LWNS = varis2D['LWNS'][:]

SEF = LHF + SHF

#clear-sky radiation variables
SWNTC = varis2D['SWNTC'][:]
SWNSC = varis2D['SWNSC'][:]
LWNTC = varis2D['LWNTC'][:]
LWNSC = varis2D['LWNSC'][:]

netSW = SWNT-SWNS
netLW = LWNS-LWNT
netSWC = SWNTC-SWNSC
netLWC = LWNSC-LWNTC

CRH = PW/SWVP
#3D variables
print 'loading 3D variables'
qv = varis3D['QV'][:]
qv = qv*1e-3
T = varis3D['TABS'][:]
w = varis3D['W'][:]

#1D variables
rhoz = p/(c.Rd*np.mean(np.mean(T, axis=3), axis=2))


#averaging time period for 2D & 3D fields
ntave2D=24
ntave3D=4

nt3D = t3D.size
nt2D = t2D.size
nz = z.size
nx = x.size
ny = y.size

dt3D = np.diff(t3D)[0]

dt3D = dt3D*(3600*24) #convert from days to seconds

#print 'calculating temporal averages'
##3D fields (6-hourly output)
##reshape to take daily averages
#qv_tmp = qv.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
#T_tmp = T.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
#w_tmp = w.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
#rho_tmp = rhoz.reshape(nt3D/ntave3D, ntave3D, nz)
#
#qv_tave = np.mean(qv_tmp, axis=1)
#T_tave = np.mean(T_tmp, axis=1)
#w_tave = np.mean(w_tmp, axis=1)
#rhoz_tave = np.mean(rho_tmp, axis=1)
#
##2D fields (hourly output)
##reshape to take daily averages
#PW_tmp = PW.reshape(nt2D/ntave2D, ntave2D, nx, ny)
#CRH_tmp = CRH.reshape(nt2D/ntave2D, ntave2D, nx, ny)
#LHF_tmp = LHF.reshape(nt2D/ntave2D, ntave2D, nx, ny)
#SHF_tmp = SHF.reshape(nt2D/ntave2D, ntave2D, nx, ny)
#IWP_tmp =  IWP.reshape(nt2D/ntave2D, ntave2D, nx, ny)
#netSW_tmp = netSW.reshape(nt2D/ntave2D, ntave2D, nx, ny)
#netLW_tmp = netLW.reshape(nt2D/ntave2D, ntave2D, nx, ny)
#netSWC_tmp = netSWC.reshape(nt2D/ntave2D, ntave2D, nx, ny)
#netLWC_tmp = netLWC.reshape(nt2D/ntave2D, ntave2D, nx, ny)
#
#PW_tave = np.mean(PW_tmp, axis=1)
#CRH_tave = np.mean(CRH_tmp, axis=1)
#LHF_tave = np.mean(LHF_tmp, axis=1)
#SHF_tave = np.mean(SHF_tmp, axis=1)
#IWP_tave = np.mean(IWP_tmp, axis=1)
#netSW_tave = np.mean(netSW_tmp, axis=1)
#netLW_tave = np.mean(netLW_tmp, axis=1)
#netSWC_tave = np.mean(netSWC_tmp, axis=1)
#netLWC_tave = np.mean(netLWC_tmp, axis=1)
#SEF_tave = LHF_tave + SHF_tave
#
#netSWcld_tave = netSW_tave - netSWC_tave
#netLWcld_tave = netLW_tave - netLWC_tave

netSWcld = netSW - netSWC
netLWcld = netLW - netLWC

z3D = np.ones((x.size, y.size, z.size))
z3D[:,:,:]=z
z3D = np.transpose(z3D)

#width of block in units of dx
#dx = 3 km, so db = 16*3 = 48 km
db=16

#times=np.arange(0,130)
times=t3D
nblocks = (nx/db)*(ny/db)

CRHs = np.zeros((times.size, nblocks))
corr_SEFs = np.zeros(CRHs.shape)
corr_netLWs = np.zeros(CRHs.shape)
corr_netSWs = np.zeros(CRHs.shape)
corr_netSWclds = np.zeros(CRHs.shape)
corr_netLWclds = np.zeros(CRHs.shape)
corr_netSWCs = np.zeros(CRHs.shape)
corr_netLWCs = np.zeros(CRHs.shape)
corr_convs = np.zeros(CRHs.shape)

hfhatprimes = np.zeros((times.size, nx, ny))
hfhatprimesort = np.zeros(CRHs.shape)



#loop over times to calculate FMSE feedback terms
for i3, ti in enumerate(times):
    
    i = (ntave2D/ntave3D) * i3
    

    
    #MSE and FMSE calculations
    print 'calculating vertically integrated MSE and FMSE'
    print ti
    
     #is it possible to get MSE and FMSE from 2D fields? 
    #problem: T is coming from 3D field, which is output 6-hourly as opposed to hourly
    #this is an inconsistency when daily averaging is performed.
    That = vertint(T[i3,:,:,:], p)
    zhat = vertint(z3D, p)
    qvhat = c.rhol*1e-3*PW[i,:,:]
    qihat = c.rhoi*1e-3*IWP[i,:,:]
    
    #vertically integrated MSE calculation
    hhat = c.cp*That + c.g*zhat + c.lv0*qvhat
    
    #zonal mean of vertical integral of MSE
    hhatbar = np.mean(hhat)
    #deviation from zonal mean of MSE
    hhatprime = hhat - hhatbar
    varhhatprime = np.var(hhatprime)
    
    #vertically integrated FMSE 
    hfhat = c.cp*That + c.g*zhat + c.lv0*qvhat - c.lf*qihat
    #zonal mean of vertical integral of FMSE
    hfhatbar = np.mean(hfhat)
    #deviation from zonal mean of FMSE
    hfhatprime = hfhat - hfhatbar
    varhfhatprime = np.var(hfhatprime)
    
    hfhatprimes[i3,:,:] = hfhatprime
    
    

    
    #calculate hhat'(SEF' + NetSW' + NetLW' )
    
    SEFprime = SEF[i,:,:] - np.mean(SEF[i,:,:])
    netLWprime = netLW[i,:,:] - np.mean(netLW[i,:,:])
    netSWprime = netSW[i,:,:] - np.mean(netSW[i,:,:])
    netSWcldprime = netSWcld[i,:,:] - np.mean(netSWcld[i,:,:])
    netLWcldprime = netLWcld[i,:,:] - np.mean(netLWcld[i,:,:])
    netSWCprime = netSWC[i,:,:] - np.mean(netSWC[i,:,:])
    netLWCprime = netLWC[i,:,:] - np.mean(netLWC[i,:,:])
    

    corr_SEF = (hfhatprime*SEFprime)
    corr_netLW = (hfhatprime*netLWprime)
    corr_netSW = (hfhatprime*netSWprime)
    corr_netSWcld = (hfhatprime*netSWcldprime)
    corr_netLWcld = (hfhatprime*netLWcldprime)
    corr_netSWC = (hfhatprime*netSWCprime)
    corr_netLWC = (hfhatprime*netLWCprime)
    
    print 'CRH block sorting'
    #CRH_SEFsort = blocksort2D(CRH[i,:,:], corr_SEF, db)
    #CRH_netLWsort = blocksort2D(CRH[i,:,:], corr_netLW, db)
    #CRH_netSWsort = blocksort2D(CRH[i,:,:], corr_netSW, db)
    #CRH_netSWcldsort = blocksort2D(CRH[i,:,:], corr_netSWcld, db)
    #CRH_netLWcldsort = blocksort2D(CRH[i,:,:], corr_netLWcld, db)
    #CRH_netSWCsort = blocksort2D(CRH[i,:,:], corr_netSWC, db)
    #CRH_netLWCsort = blocksort2D(CRH[i,:,:], corr_netLWC, db)
    #CRH_hfhatprimesort = blocksort2D(CRH[i, :, :], hfhatprime, db)
    
    CRHblock = blockave2D(CRH[i,:,:], db).flatten()
    CRH_SEFblock = blockave2D(corr_SEF, db).flatten()
    CRH_netLWblock = blockave2D(corr_netLW, db).flatten()
    CRH_netSWblock = blockave2D(corr_netSW, db).flatten()
    CRH_netSWcldblock = blockave2D(corr_netSWcld, db).flatten()
    CRH_netLWcldblock = blockave2D(corr_netLWcld, db).flatten()
    CRH_netSWCblock = blockave2D( corr_netSWC, db).flatten()
    CRH_netLWCblock = blockave2D(corr_netLWC, db).flatten()
    CRH_hfhatprimeblock = blockave2D(hfhatprime, db).flatten()
    
    
        
            
    
        

    
    CRH_sort = [y for y, x in sorted(zip(CRHblock, CRH_SEFblock))]
    CRH_SEFsort = [x for y, x in sorted(zip(CRHblock, CRH_SEFblock))]    
    
    CRH_netLWsort = [x for y, x in sorted(zip(CRHblock, CRH_netLWblock))]  
    
    CRH_netSWsort = [x for y, x in sorted(zip(CRHblock, CRH_netSWblock))]  
    
    CRH_netSWcldsort = [x for y, x in sorted(zip(CRHblock, CRH_netSWcldblock))]

    CRH_netLWcldsort = [x for y, x in sorted(zip(CRHblock, CRH_netLWcldblock))]
    
    CRH_netSWCsort = [x for y, x in sorted(zip(CRHblock, CRH_netSWCblock))] 
    
    CRH_netLWCsort = [x for y, x in sorted(zip(CRHblock, CRH_netLWCblock))] 
    
    CRH_hfhatprimesort = [x for y, x in sorted(zip(CRHblock, CRH_hfhatprimeblock))] 

  
         
        
        
    

                
    #oCRH_SEFblock = dict(zip(CRHblock, CRH_SEFblock))
    #CRH_SEFsort = collections.OrderedDict(sorted(oCRH_SEFblock.items()))
    #oCRH_netLWblock = dict(zip(CRHblock, CRH_netLWblock))
    #CRH_netLWsort = collections.OrderedDict(sorted(oCRH_netLWblock.items()))
    #oCRH_netSWblock = dict(zip(CRHblock, CRH_netSWblock))
    #CRH_netSWsort = collections.OrderedDict(sorted(oCRH_netSWblock.items()))
    #oCRH_netSWcldblock = dict(zip(CRHblock, CRH_netSWcldblock))
    #CRH_netSWcldsort = collections.OrderedDict(sorted(oCRH_netSWcldblock.items()))
    #oCRH_netLWcldblock = dict(zip(CRHblock, CRH_netLWcldblock))
    #CRH_netLWcldsort = collections.OrderedDict(sorted(oCRH_netLWcldblock.items()))   
    #oCRH_netSWCblock = dict(zip(CRHblock, CRH_netSWCblock))
    #CRH_netSWCsort = collections.OrderedDict(sorted(oCRH_netSWCblock.items()))
    #oCRH_netLWCblock = dict(zip(CRHblock, CRH_netLWCblock))
    #CRH_netLWCsort = collections.OrderedDict(sorted(oCRH_netLWCblock.items()))
    #oCRH_hfhatprimeblock = dict(zip(CRHblock, CRH_hfhatprimeblock))
    #CRH_hfhatprimesort = collections.OrderedDict(sorted(oCRH_hfhatprimeblock.items()))
    
    CRHs[i3,:] = CRH_sort
    corr_SEFs[i3,:] = CRH_SEFsort
    corr_netLWs[i3,:] = CRH_netLWsort
    corr_netSWs[i3,:] = CRH_netSWsort
    corr_netSWclds[i3,:] = CRH_netSWcldsort
    corr_netLWclds[i3,:] = CRH_netLWcldsort
    corr_netSWCs[i3,:] = CRH_netSWCsort
    corr_netLWCs[i3,:] = CRH_netLWCsort
    hfhatprimesort[i3,:] = CRH_hfhatprimesort
    
 
    #
    #CRHs[i3,:] = CRH_SEFsort.keys()
    #corr_SEFs[i3,:] = CRH_SEFsort.values()
    #corr_netLWs[i3,:] = CRH_netLWsort.values()
    #corr_netSWs[i3,:] = CRH_netSWsort.values()
    #corr_netSWclds[i3,:] = CRH_netSWcldsort.values()
    #corr_netLWclds[i3,:] = CRH_netLWcldsort.values()
    #corr_netSWCs[i3,:] = CRH_netSWCsort.values()
    #corr_netLWCs[i3,:] = CRH_netLWCsort.values()
    #hfhatprimesort[i3,:] = CRH_hfhatprimesort.values()

    #calculate correlation of h^' with horizontal convergence of FMSE
    if ti > times[0]:
        ddthfhatprime = (hfhatprimes[i3,:,:] - hfhatprimes[i3-1,:,:])/dt3D
        conv = ddthfhatprime  - (SEFprime + netLWprime + netSWprime) 
        corr_conv = (hfhatprime*conv)
        corr_convblock = blockave2D(corr_conv, db).flatten()
        #oCRH_corr_convblock = dict(zip(CRHblock, corr_convblock))
        CRH_corr_convsort = [x for y, x in sorted(zip(CRHblock, corr_convblock))] 
        #CRH_convsort = blocksort2D(CRH[i,:,:], corr_conv, db)
        corr_convs[i3,:] = CRH_corr_convsort

#going to want to contour correlations in 2D (CRH rank, time)
CRHranks = np.arange(np.size(CRHs[0,:]))
rankss, tt = np.meshgrid(CRHranks, times)
#convert tt to days
#tt = tt+t3D[-1]

####### PLOTTING #######

#plot MSE correlation with different MSE budget terms in 2D (CRH rank, time)
plt.figure(1)
f, axarr = plt.subplots(7,1)
vmin = -1
vmax = 1

hfhatprimesortvar = np.multiply(hfhatprimesort, hfhatprimesort)

hfhatvarbar_t = np.mean(hfhatprimesortvar, axis=1)

hfhatvarbar = np.zeros(hfhatprimesortvar.shape).T
hfhatvarbar[:,:] = hfhatvarbar
hfhatvarbar = hfhatvarbar.T

levels = np.linspace(-1,1,100)

cf = axarr[0].contourf(rankss, tt, (corr_SEFs/hfhatvarbar)*(3600*24), levels, cmap=cm.RdBu_r)
axarr[0].set_title(r'$\hat h^{\prime} \mathrm{SEF}^{\prime}$')
cf = axarr[1].contourf(rankss, tt, (corr_netLWs/hfhatvarbar)*(3600*24), levels, cmap=cm.RdBu_r)
axarr[1].set_title(r'$\hat h^{\prime} \mathrm{netLW}^{\prime}$')
cf = axarr[2].contourf(rankss, tt, (corr_netSWs/hfhatvarbar)*(3600*24), levels, cmap=cm.RdBu_r)
axarr[2].set_title(r'$\hat h^{\prime} \mathrm{netSW}^{\prime}$')
cf = axarr[3].contourf(rankss, tt, (corr_netSWclds/hfhatvarbar)*(3600*24), levels, cmap=cm.RdBu_r)
axarr[3].set_title(r'$\hat h^{\prime} \mathrm{netSW}_{cld}^{\prime}$')
cf = axarr[4].contourf(rankss, tt, (corr_netLWclds/hfhatvarbar)*(3600*24), levels, cmap=cm.RdBu_r)
axarr[4].set_title(r'$\hat h^{\prime} \mathrm{netLW}_{cld}^{\prime}$')
cf = axarr[5].contourf(rankss, tt, (corr_netSWCs/hfhatvarbar)*(3600*24), levels, cmap=cm.RdBu_r)
axarr[5].set_title(r'$\hat h^{\prime} \mathrm{netLW}_{clear}^{\prime}$')
cf = axarr[6].contourf(rankss, tt, (corr_netLWCs/hfhatvarbar)*(3600*24), levels, cmap=cm.RdBu_r)
axarr[6].set_title(r'$\hat h^{\prime} \mathrm{netSW}_{clear}^{\prime}$')
cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
f.text(0.47, 0.06, 'CRH rank', ha='center', va='center')
f.text(0.08, 0.5, 'time (days)', ha='center', va='center', rotation='vertical')
f.subplots_adjust(right=0.8)
cb = f.colorbar(cf, cax=cbar_ax)
#TODO: add superior x,y axis labels
#plt.ylabel('days')
cb.set_label('FMSE anomaly amplification (per day)')
plt.savefig(fout + 'CRHsort_feedbacks_day{0}to{1}.pdf'.format(tt[0, 0], tt[-1,0]))
plt.close()

plt.figure(2)
vmin =-1
vmax = 1
cf = plt.contourf(rankss, tt, (corr_convs/hfhatvarbar)*(3600*24), levels, cmap=cm.RdBu_r, vmin=vmin, vmax=vmax)
plt.title(r'$-\hat h^{\prime} \nabla _H  \cdot{\hat{\vec{u}}h}$')
plt.xlabel('CRH rank')
plt.ylabel('time (days)')
#plt.ylim(-1,1)
plt.colorbar(label='FMSE anomaly amplification (per day)')
plt.savefig(fout + 'CRHsort_convfeedback_day{0}to{1}.pdf'.format(tt[0, 0], tt[-1,0]))
plt.close()

#plot domain mean correlations

corr_SEFbar = np.mean(corr_SEFs, axis=1)
corr_netLWbar = np.mean(corr_netLWs, axis=1)
corr_netSWbar = np.mean(corr_netSWs, axis=1)
corr_netSWcldbar = np.mean(corr_netSWclds, axis=1)
corr_netLWcldbar = np.mean(corr_netLWclds, axis=1)
corr_netSWCbar = np.mean(corr_netSWCs, axis=1)
corr_netLWCbar = np.mean(corr_netLWCs, axis=1)
corr_convbar = np.mean(corr_convs, axis=1)

nave3D=4
w = nave3D*5
conv_smooth = moving_average((corr_convbar/hfhatvarbar), w)

plt.figure(2)
plt.plot(tt[:,0], (corr_SEFbar/hfhatvarbar)*(3600*24), label=r'$\overline{\hat h^{\prime} \mathrm{SEF}^{\prime}}$')
plt.plot(tt[:,0], (corr_netLWbar/hfhatvarbar)*(3600*24), label=r'$\overline{\hat h^{\prime} \mathrm{netLW}^{\prime}}$')
plt.plot(tt[:,0], (corr_netSWbar/hfhatvarbar)*(3600*24), label=r'$\overline{\hat h^{\prime} \mathrm{netSW}^{\prime}}$')
#plt.plot(tt[:,0], corr_netSWcldbar*(3600*24), label=r'$\overline{\hat h^{\prime} \mathrm{netSW}_{cld}^{\prime}}$')
#plt.plot(tt[:,0], corr_netLWcldbar*(3600*24), label=r'$\overline{\hat h^{\prime} \mathrm{netLW}_{cld}^{\prime}}$')
#plt.plot(tt[:,0], corr_netLWCbar*(3600*24), label=r'$\overline{\hat h^{\prime} \mathrm{netLW}_{clear}^{\prime}}$')
#plt.plot(tt[:,0], corr_netSWCbar*(3600*24), label=r'$\overline{\hat h^{\prime} \mathrm{netSW}_{clear}^{\prime}}$')
plt.plot(tt[:,0], (corr_convbar/hfhatvarbar)*(3600*24), 'g--', label=r'$- \overline{\hat h^{\prime} \nabla _H \cdot{\hat{\vec{u}}h}}$')
plt.plot(tt[:,w/2:-w/2+1], conv_smooth*(3600*24), 'g', label=r'smoothed convergence')
plt.title(r'FMSE variance budget terms, normalized by domain mean FMSE ')
plt.xlabel('time (days)')
plt.ylabel('rate of change of FMSE variance (per day)')
plt.ylim(-1,1)
plt.legend(loc='best')
plt.savefig(fout + 'mean_feedbacks_day{0}to{1}.pdf'.format(tt[0, 0], tt[-1,0]))
plt.close()

plt.figure(3)
plt.plot(tt[:,0], (corr_netSWcldbar/hfhatvarbar)*(3600*24), label=r'$\overline{\hat h^{\prime} \mathrm{netSW}_{cld}^{\prime}}$')
plt.plot(tt[:,0], (corr_netLWcldbar/hfhatvarbar)*(3600*24), label=r'$\overline{\hat h^{\prime} \mathrm{netLW}_{cld}^{\prime}}$')
plt.plot(tt[:,0], (corr_netLWCbar/hfhatvarbar)*(3600*24), label=r'$\overline{\hat h^{\prime} \mathrm{netLW}_{clear}^{\prime}}$')
plt.plot(tt[:,0], (corr_netSWCbar/hfhatvarbar)*(3600*24), label=r'$\overline{\hat h^{\prime} \mathrm{netSW}_{clear}^{\prime}}$')
plt.title(r'FMSE variance budget terms (radiation), normalized by domain mean FMSE')
plt.xlabel('time (days)')
plt.ylabel('rate of change of FMSE variance (per day)')
plt.legend(loc='best')
plt.ylim(-1,1)
plt.savefig(fout + 'mean_radfeedbacks_day{0}to{1}.pdf'.format(tt[0, 0], tt[-1,0]))
plt.close()

#plot the evolution of the variance of FMSE, sorted by CRH

plt.figure(4)
plt.contourf(rankss, tt, hfhatprimesortvar, 60, cmap=cm.RdBu_r)
plt.colorbar()
plt.contour(rankss, tt, hfhatprimesort, levels=[0], colors='k')
plt.ylabel('time (days)')
plt.xlabel('CRH rank')
plt.title(r'$\hat h^{\prime 2} (J^2/m^4)$')
plt.savefig(fout + 'CRHsort_FMSEvar_day{0}to{1}.pdf'.format(tt[0, 0], tt[-1,0]))
plt.close()














