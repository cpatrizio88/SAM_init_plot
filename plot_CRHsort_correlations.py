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

c = constants()

p_s = 1015
T_s = 302

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (18, 12)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
#fout = '/Users/cpatrizio/figures/SST302/SAM_aggr130days_768km_64vert_ubarzero_CRHSORT/'
fout = '/Users/cpatrizio/figures/SST302/1536km_SAM_aggr140days_64vert_ubarzero_CRHSORT'


#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day090to140*302K.nc')[0]
nc_in2D = glob.glob(fpath2D + '*512x512*3000m*day090to140*302K.nc')[0]

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
PW = varis2D['PW'][:]
IWP = varis2D['IWP'][:]
LHF = varis2D['LHF'][:]
SHF = varis2D['SHF'][:]
SWVP = varis2D['SWVP'][:]
SWNT = varis2D['SWNT'][:]
SWNS = varis2D['SWNS'][:]
LWNT = varis2D['LWNT'][:]
LWNS = varis2D['LWNS'][:]

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

print 'calculating temporal averages'
#3D fields (6-hourly output)
#reshape to take daily averages
qv_tmp = qv.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
T_tmp = T.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
w_tmp = w.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
rho_tmp = rhoz.reshape(nt3D/ntave3D, ntave3D, nz)

qv_tave = np.mean(qv_tmp, axis=1)
T_tave = np.mean(T_tmp, axis=1)
w_tave = np.mean(w_tmp, axis=1)
rhoz_tave = np.mean(rho_tmp, axis=1)

#2D fields (hourly output)
#reshape to take daily averages
PW_tmp = PW.reshape(nt2D/ntave2D, ntave2D, nx, ny)
CRH_tmp = CRH.reshape(nt2D/ntave2D, ntave2D, nx, ny)
LHF_tmp = LHF.reshape(nt2D/ntave2D, ntave2D, nx, ny)
SHF_tmp = SHF.reshape(nt2D/ntave2D, ntave2D, nx, ny)
IWP_tmp =  IWP.reshape(nt2D/ntave2D, ntave2D, nx, ny)
netSW_tmp = netSW.reshape(nt2D/ntave2D, ntave2D, nx, ny)
netLW_tmp = netLW.reshape(nt2D/ntave2D, ntave2D, nx, ny)
netSWC_tmp = netSWC.reshape(nt2D/ntave2D, ntave2D, nx, ny)
netLWC_tmp = netLWC.reshape(nt2D/ntave2D, ntave2D, nx, ny)

PW_tave = np.mean(PW_tmp, axis=1)
CRH_tave = np.mean(CRH_tmp, axis=1)
LHF_tave = np.mean(LHF_tmp, axis=1)
SHF_tave = np.mean(SHF_tmp, axis=1)
IWP_tave = np.mean(IWP_tmp, axis=1)
netSW_tave = np.mean(netSW_tmp, axis=1)
netLW_tave = np.mean(netLW_tmp, axis=1)
netSWC_tave = np.mean(netSWC_tmp, axis=1)
netLWC_tave = np.mean(netLWC_tmp, axis=1)
SEF_tave = LHF_tave + SHF_tave

netSWcld_tave = netSW_tave - netSWC_tave
netLWcld_tave = netLW_tave - netLWC_tave

z3D = np.ones((x.size, y.size, z.size))
z3D[:,:,:]=z
z3D = np.transpose(z3D)

#width of block in units of dx
#dx = 3 km, so db = 16*3 = 48 km
db=16

times=np.arange(0,130)
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
for i, ti in enumerate(times):
    
    #MSE and FMSE calculations
    print 'calculating vertically integrated MSE and FMSE'
    print ti
    
     #is it possible to get MSE and FMSE from 2D fields? 
    #problem: T is coming from 3D field, which is output 6-hourly as opposed to hourly
    #this is an inconsistency when daily averaging is performed.
    That = vertint(T_tave[ti,:,:,:], p)
    zhat = vertint(z3D, p)
    qvhat = c.rhol*1e-3*PW_tave[ti,:,:]
    qihat = c.rhoi*1e-3*IWP_tave[ti,:,:]
    
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
    
    hfhatprimes[ti,:,:] = hfhatprime
    

    
    #calculate hhat'(SEF' + NetSW' + NetLW' )
    
    SEFprime = SEF_tave[ti,:,:] - np.mean(SEF_tave[ti,:,:])
    netLWprime = netLW_tave[ti,:,:] - np.mean(netLW_tave[ti,:,:])
    netSWprime = netSW_tave[ti,:,:] - np.mean(netSW_tave[ti,:,:])
    netSWcld = netSWcld_tave[ti,:,:] - np.mean(netSWcld_tave[ti,:,:])
    netLWcld = netLWcld_tave[ti,:,:] - np.mean(netLWcld_tave[ti,:,:])
    netSWC = netSWC_tave[ti,:,:] - np.mean(netSWC_tave[ti,:,:])
    netLWC = netLWC_tave[ti,:,:] - np.mean(netLWC_tave[ti,:,:])
    

    corr_SEF = (hfhatprime*SEFprime)/varhfhatprime
    corr_netLW = (hfhatprime*netLWprime)/varhfhatprime
    corr_netSW = (hfhatprime*netSWprime)/varhfhatprime
    corr_netSWcld = (hfhatprime*netSWcld)/varhfhatprime
    corr_netLWcld = (hfhatprime*netLWcld)/varhfhatprime
    corr_netSWC = (hfhatprime*netSWC)/varhfhatprime
    corr_netLWC = (hfhatprime*netLWC)/varhfhatprime
    
    print 'CRH block sorting'
    CRH_SEFsort = blocksort2D(CRH_tave[ti,:,:], corr_SEF, db)
    CRH_netLWsort = blocksort2D(CRH_tave[ti,:,:], corr_netLW, db)
    CRH_netSWsort = blocksort2D(CRH_tave[ti,:,:], corr_netSW, db)
    CRH_netSWcldsort = blocksort2D(CRH_tave[ti,:,:], corr_netSWcld, db)
    CRH_netLWcldsort = blocksort2D(CRH_tave[ti,:,:], corr_netLWcld, db)
    CRH_netSWCsort = blocksort2D(CRH_tave[ti,:,:], corr_netSWC, db)
    CRH_netLWCsort = blocksort2D(CRH_tave[ti,:,:], corr_netLWC, db)
    CRH_hfhatprimesort = blocksort2D(CRH_tave[ti, :, :], hfhatprime, db)
 
    
    CRHs[i,:] = CRH_SEFsort.keys()
    corr_SEFs[i,:] = CRH_SEFsort.values()
    corr_netLWs[i,:] = CRH_netLWsort.values()
    corr_netSWs[i,:] = CRH_netSWsort.values()
    corr_netSWclds[i,:] = CRH_netSWcldsort.values()
    corr_netLWclds[i,:] = CRH_netLWcldsort.values()
    corr_netSWCs[i,:] = CRH_netSWCsort.values()
    corr_netLWCs[i,:] = CRH_netLWCsort.values()
    hfhatprimesort[i,:] = CRH_hfhatprimesort.values()

    #calculate correlation of h^' with horizontal convergence of FMSE
    if ti > times[0]:
        ddthfhatprime = (hfhatprimes[ti,:,:] - hfhatprimes[ti-1,:,:])/dt3D
        conv = ddthfhatprime  - (SEFprime + netLWprime + netSWprime) 
        corr_conv = (hfhatprime*conv)/varhfhatprime
        CRH_convsort = blocksort2D(CRH_tave[ti,:,:], corr_conv, db)
        corr_convs[i,:] = CRH_convsort.values()

#going to want to contour correlations in 2D (CRH rank, time)
CRHranks = np.arange(np.size(CRHs[0,:]))
rankss, tt = np.meshgrid(CRHranks, times)
#convert tt to days
#tt = tt+t3D[-1]

####### PLOTTING #######

#plot MSE correlation with different MSE budget terms in 2D (CRH rank, time)
plt.figure(1)
f, axarr = plt.subplots(7,1)
vmin = -1e-5
vmax = 1e-5
cf = axarr[0].contourf(rankss, tt, corr_SEFs, 60, cmap=cm.RdBu_r, vmin=vmin, vmax=vmax)
axarr[0].set_title(r'$\overline{\hat h^{\prime} \mathrm{SEF}^{\prime}}$')
cf = axarr[1].contourf(rankss, tt, corr_netLWs, 60, cmap=cm.RdBu_r, vmin=vmin, vmax=vmax)
axarr[1].set_title(r'$\overline{\hat h^{\prime} \mathrm{netLW}^{\prime}}$')
cf = axarr[2].contourf(rankss, tt, corr_netSWs, 60, cmap=cm.RdBu_r, vmin=vmin, vmax=vmax)
axarr[2].set_title(r'$\overline{\hat h^{\prime} \mathrm{netSW}^{\prime}}$')
cf = axarr[3].contourf(rankss, tt, corr_netSWclds, 60, cmap=cm.RdBu_r, vmin=vmin, vmax=vmax)
axarr[3].set_title(r'$\overline{\hat h^{\prime} \mathrm{netSW}_{cld}^{\prime}}$')
cf = axarr[4].contourf(rankss, tt, corr_netLWclds, 60, cmap=cm.RdBu_r, vmin=vmin, vmax=vmax)
axarr[4].set_title(r'$\overline{\hat h^{\prime} \mathrm{netLW}_{cld}^{\prime}}$')
cf = axarr[5].contourf(rankss, tt, corr_netSWCs, 60, cmap=cm.RdBu_r, vmin=vmin, vmax=vmax)
axarr[5].set_title(r'$\overline{\hat h^{\prime} \mathrm{netLW}_{clear}^{\prime}}$')
cf = axarr[6].contourf(rankss, tt, corr_netLWCs, 60, cmap=cm.RdBu_r, vmin=vmin, vmax=vmax)
axarr[6].set_title(r'$\overline{\hat h^{\prime} \mathrm{netSW}_{clear}^{\prime}}$')
cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
f.text(0.47, 0.06, 'CRH rank', ha='center', va='center')
f.text(0.08, 0.5, 'time (days)', ha='center', va='center', rotation='vertical')
f.subplots_adjust(right=0.8)
cb = f.colorbar(cf, cax=cbar_ax)
#TODO: add superior x,y axis labels
#plt.ylabel('days')
cb.set_label('correlation')
plt.savefig(fout + 'CRHsort_feedbacks_day{0}to{1}.pdf'.format(tt[0, 0], tt[-1,0]))
plt.close()

plt.figure(2)
vmin = -2e-4
vmax = 2e-4
cf = plt.contourf(rankss, tt, corr_convs, 60, cmap=cm.RdBu_r, vmin=vmin, vmax=vmax)
plt.title(r'$- \overline{\hat h^{\prime} \nabla _H  \cdot{\hat{\vec{u}}h}}$')
plt.xlabel('CRH rank')
plt.ylabel('time (days)')
plt.colorbar()
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

plt.figure(2)
plt.plot(tt[:,0], corr_SEFbar, label=r'$\overline{\hat h^{\prime} \mathrm{SEF}^{\prime}}$')
plt.plot(tt[:,0], corr_netLWbar, label=r'$\overline{\hat h^{\prime} \mathrm{netLW}^{\prime}}$')
plt.plot(tt[:,0], corr_netSWbar, label=r'$\overline{\hat h^{\prime} \mathrm{netSW}^{\prime}}$')
plt.plot(tt[:,0], corr_netSWcldbar, label=r'$\overline{\hat h^{\prime} \mathrm{netSW}_{cld}^{\prime}}$')
plt.plot(tt[:,0], corr_netLWcldbar, label=r'$\overline{\hat h^{\prime} \mathrm{netLW}_{cld}^{\prime}}$')
plt.plot(tt[:,0], corr_netLWCbar, label=r'$\overline{\hat h^{\prime} \mathrm{netLW}_{clear}^{\prime}}$')
plt.plot(tt[:,0], corr_netSWCbar, label=r'$\overline{\hat h^{\prime} \mathrm{netSW}_{clear}^{\prime}}$')
plt.plot(tt[:,0], corr_convbar, 'g--', label=r'$- \overline{\hat h^{\prime} \nabla _H \cdot{\hat{\vec{u}}h}}$')
plt.title(r'correlation between vertically integrated MSE anomaly, $\hat h^{\prime}$, and budget terms')
plt.xlabel('days')
plt.ylabel('correlation')
plt.legend(loc='best')
plt.savefig(fout + 'mean_feedbacks_day{0}to{1}.pdf'.format(tt[0, 0], tt[-1,0]))
plt.close()

#plot the evolution of the variance of FMSE, sorted by CRH
hfhatprimesortvar = np.multiply(hfhatprimesort, hfhatprimesort)
plt.figure(3)
plt.contourf(rankss, tt, hfhatprimesortvar, 60, cmap=cm.RdBu_r)
plt.colorbar()
plt.contour(rankss, tt, hfhatprimesort, levels=[0], colors='k')
plt.ylabel('time (days)')
plt.xlabel('CRH rank')
plt.title(r'$\hat h^{\prime 2} (J^2/m^4)$')
plt.savefig(fout + 'CRHsort_FMSEvar_day{0}to{1}.pdf'.format(tt[0, 0], tt[-1,0]))
plt.close()














