from netCDF4 import Dataset
import site
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
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/figures/SST302/SAM_aggrday90to130_768km_64vert_ubarzero_3D/'

print 'loading 3D variables'
nc_in3D = glob.glob(fpath3D + '*256x256*3000m*day90to130*302K.nc')[0]
print 'loading 2D variables'
nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

nc_data3D = Dataset(nc_in3D)
nc_data2D = Dataset(nc_in2D)
varis3D = nc_data3D.variables
varis2D = nc_data2D.variables

t3D = varis3D['time'][:]
t2D = varis2D['time'][:]
x = varis3D['x'][:]
y = varis3D['y'][:]
xx, yy = np.meshgrid(x, y)
z = varis3D['z'][:]
p = varis3D['p'][:]
p = p*1e2
#2D variables
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
qv = varis3D['QV'][:]
qv = qv*1e-3
T = varis3D['TABS'][:]

ntave2D=24
ntave3D=4

nt3D = t3D.size
nt2D = t2D.size
nz = z.size
nx = x.size
ny = y.size

print 'calculating temporal averages'
#3D fields (6-hourly output)
#reshape to take daily averages
qv_tmp = qv.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)
T_tmp = T.reshape(nt3D/ntave3D, ntave3D, nz, nx, ny)

qv_tave = np.mean(qv_tmp, axis=1)
T_tave = np.mean(T_tmp, axis=1)

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

#LOOP OVER DAYS LATER TO SEE TEMPORAL EVOLUTION OF STEADY STATE ?
times=np.arange(-40,0)
nblocks = (nx/db)*(ny/db)

CRHs = np.zeros((times.size, nblocks))
corr_SEFs = np.zeros(CRHs.shape)
corr_netLWs = np.zeros(CRHs.shape)
corr_netSWs = np.zeros(CRHs.shape)
corr_netSWclds = np.zeros(CRHs.shape)
corr_netLWclds = np.zeros(CRHs.shape)

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
    
    #calculate hhat'(SEF' + NetSW' + NetLW' )
    
    SEFprime = SEF_tave[ti,:,:] - np.mean(SEF_tave[ti,:,:])
    netLWprime = netLW_tave[ti,:,:] - np.mean(netLW_tave[ti,:,:])
    netSWprime = netSW_tave[ti,:,:] - np.mean(netSW_tave[ti,:,:])
    netSWcld = netSWcld_tave[ti,:,:] - np.mean(netSWcld_tave[ti,:,:])
    netLWcld = netLWcld_tave[ti,:,:] - np.mean(netLWcld_tave[ti,:,:])
    
    
    corr_SEF = (hfhatprime*SEFprime)/varhfhatprime
    corr_netLW = (hfhatprime*netLWprime)/varhfhatprime
    corr_netSW = (hfhatprime*netSWprime)/varhfhatprime
    corr_netSWcld = (hfhatprime*netSWcld)/varhfhatprime
    corr_netLWcld = (hfhatprime*netLWcld)/varhfhatprime
    
    print 'CRH block sorting'
    CRH_SEFsort = blocksort2D(CRH_tave[ti,:,:], corr_SEF, db)
    CRH_netLWsort = blocksort2D(CRH_tave[ti,:,:], corr_netLW, db)
    CRH_netSWsort = blocksort2D(CRH_tave[ti,:,:], corr_netSW, db)
    CRH_netSWcldsort = blocksort2D(CRH_tave[ti,:,:], corr_netSWcld, db)
    CRH_netLWcldsort = blocksort2D(CRH_tave[ti,:,:], corr_netLWcld, db)
    
    CRHs[i,:] = CRH_SEFsort.keys()
    corr_SEFs[i,:] = CRH_SEFsort.values()
    corr_netLWs[i,:] = CRH_netLWsort.values()
    corr_netSWs[i,:] = CRH_netSWsort.values()
    corr_netSWclds[i,:] = CRH_netSWcldsort.values()
    corr_netLWclds[i,:] = CRH_netLWcldsort.values()

CRHranks = np.arange(np.size(CRHs[0,:]))

rankss, tt = np.meshgrid(CRHranks, times)
tt = tt+t2D[-1]

#test blocksort3D. working properly - do temporal averges?
CRH_qvsort = blocksort3D(CRH_tave[-1,:,:], qv_tave[-1,:,:,:], db)
qvsort = CRH_qvsort[1]

CRHs_t = CRHs[-1,:]

CRHq1 = np.percentile(CRHs_t, 25)
CRHq2 = np.percentile(CRHs_t, 50)
CRHq3 = np.percentile(CRHs_t, 75)
CRHq4 = np.percentile(CRHs_t, 100)

qvq1 = qvsort[:, CRHs_t < CRHq1]
qvq2 = qvsort[:, np.bitwise_and(CRHs_t > CRHq1, CRHs_t < CRHq2)]
qvq3 = qvsort[:, np.bitwise_and(CRHs_t > CRHq2, CRHs_t < CRHq3)]
qvq4 = qvsort[:, CRHs_t > CRHq3]

qvq1bar = np.mean(qvq1, axis=1)
qvq2bar = np.mean(qvq2, axis=1)
qvq3bar = np.mean(qvq3, axis=1)
qvq4bar = np.mean(qvq4, axis=1)

#test make PW 2D (CRH, t) frequency map - working properly.
timeslong = np.arange(0, 130)
CRHslong =  np.zeros((timeslong.size, nblocks))
for i, t in enumerate(timeslong):
    print i
    CRHslong[i,:] = np.sort(blockave2D(CRH_tave[i,:,:], db).flatten())
    
ranksslong, ttlong = np.meshgrid(CRHranks, timeslong)
tcoords = ttlong.flatten()
CRHcoords = CRHslong.flatten()
heatmap, xedges, yedges = np.histogram2d(CRHcoords, np.transpose(tcoords), bins=(256,130))

extent = [CRHranks[0], CRHranks[-1], yedges[-1], yedges[0]]

plt.figure(3)
plt.imshow(np.transpose(heatmap), extent=extent, cmap=cm.RdYlBu_r)
plt.colorbar()
plt.show()




plt.figure(1)
f, axarr = plt.subplots(5,1)
cf = axarr[0].contourf(rankss, tt, corr_SEFs, 60, cmap=cm.RdBu_r)
axarr[0].set_title('corr_SEF')
cf = axarr[1].contourf(rankss, tt, corr_netLWs, 60, cmap=cm.RdBu_r)
axarr[1].set_title('corr_netLW')
cf = axarr[2].contourf(rankss, tt, corr_netSWs, 60, cmap=cm.RdBu_r)
axarr[2].set_title('corr_netSW')
cf = axarr[3].contourf(rankss, tt, corr_netSWclds, 60, cmap=cm.RdBu_r)
axarr[3].set_title('corr_netSW due to clouds')
cf = axarr[4].contourf(rankss, tt, corr_netLWclds, 60, cmap=cm.RdBu_r)
axarr[4].set_title('corr_netLW due to clouds')
cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
f.subplots_adjust(right=0.8)
f.colorbar(cf, cax=cbar_ax)
plt.show()

plt.figure(2)
plt.plot(qvq1bar, z/1e3, label='q1')
plt.plot(qvq2bar, z/1e3, label='q2')
plt.plot(qvq3bar, z/1e3, label='q3')
plt.plot(qvq4bar, z/1e3, label='q4')
plt.legend()
plt.title('average water vapor profiles sorted by CRH quartiles, at day {0}'.format(t3D[-1]))
plt.show()

#plt.contour(xx/1e3, yy/1e3, hhatprime/varhhatprime, 20, colors='k', alpha=0.5)
#plt.contourf(xx/1e3, yy/1e3, hhatprime/varhhatprime, 20, cmap=cm.RdYlBu_r, zorder=0)
#plt.title('h^''/var(h^'')')
#plt.colorbar()
#plt.show()
#plt.figure(2)
#plt.contour(xx/1e3, yy/1e3, hfhatprime/varhfhatprime, 20, colors='k', alpha=0.5)
#plt.contourf(xx/1e3, yy/1e3, hfhatprime/varhfhatprime, 20, cmap=cm.RdYlBu_r, zorder=0)
#plt.title('hf^''/var(hf^'')')
#plt.colorbar()
#plt.show()
#plt.figure(4)
#plt.plot(CRHs, corr_SEFs, label='corr_SEF')
#plt.plot(CRHs, corr_netLWs, label='corr_netLW')
#plt.plot(CRHs, corr_netSWs, label='corr_netSW')
#plt.legend()
#plt.title('correlation between total surface heat flux and column relative humidity at day={:2.1f}'.format(t2D[ti]))
#plt.show()
#plt.figure(3)
#plt.contour(xx[::db,::db]/1e3, yy[::db,::db]/1e3, CRH_blocked, 20, colors='k', alpha=0.5)
#plt.contourf(xx[::db,::db]/1e3, yy[::db,::db]/1e3, CRH_blocked, 20, cmap=cm.RdYlBu_r, zorder=0)
#plt.title('blocked')
#plt.show()





