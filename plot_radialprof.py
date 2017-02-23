from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
import SAM_init_plot.block_fns
import SAM_init_plot.misc_fns
from SAM_init_plot.misc_fns import raddist, radprof
from SAM_init_plot.block_fns import blockave2D, blockxysort2D, xysort
from scipy.optimize import curve_fit

plt.close('all')

def expfunc(x, a, b, c):
    return a*np.exp(-b*x) + c

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_RADIAL/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to140_64vert_ubarzero_RADIAL/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to150_64vert_ubarzero_RADIAL/'

#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

print 'loading .nc file'

nc_in2D = glob.glob(fpath + '*256x256*3000m*day230to250*302K.nc')[0]
#nc_in2D = glob.glob(fpath+ '*512x512*3000m*day90to140*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*1024x1024*3000m*day140to150*302K.nc')[0]



domsize=768
#domsize=1536

print 'loading 2D variables'
nc_data2D = Dataset(nc_in2D)
varis2D = nc_data2D.variables

nave=5
ntave2D=24
t=-1

aveperiod = nave*ntave2D

PW = varis2D['PW'][t-aveperiod:t,:,:]
USFC = varis2D['USFC'][t-aveperiod:,:,:]
VSFC = varis2D['VSFC'][t-aveperiod:,:,:]

speedSFC= np.sqrt(np.power(USFC, 2) + np.power(VSFC,2))
LWNT = varis2D['LWNT'][t-aveperiod:t,:,:]
LWNS = varis2D['LWNS'][t-aveperiod:t,:,:]
SWNT = varis2D['SWNT'][t-aveperiod:t,:,:]
SWNS = varis2D['SWNS'][t-aveperiod:t,:,:]
    
QNTOA = SWNT - LWNT
QNS = SWNS - LWNS

SWN = SWNT - SWNS

QN = QNTOA - QNS
    
#field = np.sqrt(np.power(USFC, 2) + np.power(VSFC, 2))
t2D = varis2D['time'][t-aveperiod:t]
x = varis2D['x'][:]
y = varis2D['y'][:]
nt2D = t2D.size
nx = x.size
ny = y.size
varnames = ['PW', 'USFC', 'LHF', 'CWP', 'Prec', 'W500', 'TB', 'ZC']
#varnames = ['W500']
#varnames = ['USFC', 'W500', 'Prec', 'LHF']
#varnames = ['PW', 'USFC', 'LHF']



#ntrunc = PW.shape[0]%ntave2D

#truncate the first few time steps to get an even number of days
#PW = PW[ntrunc:,:,:]
#QNTOA = QNTOA[ntrunc:,:,:]
#QN = QN[ntrunc:,:,:]
    
#PW_tmp = PW.reshape(nt2D/ntave2D, ntave2D, nx, ny)
#PW_tave = np.mean(PW_tmp, axis=1)

PW_tave = np.mean(PW, axis=0)

#QNTOA_tmp = QNTOA.reshape(nt2D/ntave2D, ntave2D, nx, ny)
#QNTOA_tave = np.mean(QNTOA_tmp, axis=1)

QNTOA_tave = np.mean(QNTOA, axis=0)

#QN_tmp = QN.reshape(nt2D/ntave2D, ntave2D, nx, ny)
#QN_tave = np.mean(QN_tmp, axis=1)

QN_tave = np.mean(QN, axis=0)

SWN_tave = np.mean(SWN, axis=0)

xx, yy = np.meshgrid(x, y)
times = np.arange(t2D[0], t2D[-1])
db=1
#AVERAGING PERIOD IN DAYS
nave=5

#PW_t = PW_tave[t,:,:]
#PW_t = np.mean(PW_tave[t-nave:t,:,:], axis=0)
#QNTOA_t = QNTOA_tave[t,:,:]
#QNTOA_t = np.mean(QNTOA_tave[t-nave:t,:,:], axis=0)
#QN_t = np.mean(QN_tave[t-nave:t,:,:], axis=0)

PW_blocked = blockave2D(PW_tave, db)
QNTOA_blocked=blockave2D(QNTOA_tave, db)
SWN_blocked = blockave2D(SWN_tave, db)
QN_blocked=blockave2D(QN_tave, db)

PWxy_sorted = blockxysort2D(PW_tave, xx, yy, db)
#fieldxy_sorted = blockxysort2D(field_t, xx, yy, db)

PWsort = PWxy_sorted.keys()
PWxycoords = PWxy_sorted.values()
#fieldsort = fieldxy_sorted.keys()
#fieldxycoords = fieldxy_sorted.values()

mcenter = PWxycoords[-1]

d = raddist(xx, yy, mcenter)

#put sum of PW in bins of r
#nr is r bins (indices give the value of r)
#PWrsums is sum of PW in bins of r

#EDIT: binwidth is the number of r bins.. should be consistent length with increasing
#domain size. 100 bins for 768 km domain mean 
binwidth=5e3
bins = np.arange(0, domsize*1e3, binwidth)
PWrbins, PWmeans = radprof(PW_tave, xx, yy, mcenter, binwidth)
PWrbin_centers = (PWrbins[1:] + PWrbins[:-1])/2.
#PWrprof = PWrsums/PWnr
#PWrs = np.arange(0, len(PWnr))

QNTOArbins, QNTOAmeans = radprof(QNTOA_tave, xx, yy, mcenter, binwidth)
QNTOArbin_centers = (QNTOArbins[1:] + QNTOArbins[:-1])/2.
#QNTOArprof = QNTOAsums/QNTOAnr
#QNTOArs = np.arange(0, len(QNTOAnr))

QNrbins, QNmeans = radprof(QN_tave, xx, yy, mcenter, binwidth)
QNrbin_centers = (QNrbins[1:] + QNrbins[:-1])/2.

SWNrbins, SWNmeans = radprof(SWN_tave, xx, yy, mcenter, binwidth)
SWNrbin_centers = (SWNrbins[1:] + SWNrbins[:-1])/2.
 
a=1.5 
moist_edgePW = np.mean(PW_tave) + a*np.std(PW_tave)  

for varname in varnames:
    print varname
        
    vari = varis2D[varname]
    field = varis2D[varname][t-aveperiod:t,:,:]
     #calculate speed instead of zonal wind velocity at surface
    if varname == 'USFC': 
        field=speedSFC
    #field = field[ntrunc:,:,:]
    #field_tmp = field.reshape(nt2D/ntave2D, ntave2D, nx, ny)
    #field_tave = np.mean(field_tmp, axis=1)
    field_tave = np.mean(field, axis=0)
    
    #field_t = field_tave[t,:,:]
    #average over nave days 
    #field_t = np.mean(field_tave[t-nave:t,:,:], axis=0)
    
    field_blocked = blockave2D(field_tave, db)
    
    #fieldnr, fieldsums = radprof(field_t, xx, yy, mcenter)
    fieldrbins, fieldmeans = radprof(field_tave, xx, yy, mcenter, bins)
    fieldrbin_centers = (fieldrbins[1:] + fieldrbins[:-1])/2.
    #fieldprof = fieldsums/fieldnr
    
    #fieldrs = np.arange(0, len(fieldnr))
    
    if varname == 'CWP':
        cvals = np.arange(0, 1.1, .02)
    elif varname == 'ZC':
        cvals = np.arange(0, 13, 0.25)
    elif varname == 'Prec':
        cvals = np.arange(0, 250, 5)
    elif varname == 'W500':
        cvals = np.arange(0, 0.45, 0.025)
    else:
        cvals = 60
    
    
    fig = plt.figure(1)
    ax = fig.gca()
    plt.contourf(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, field_blocked, cvals, cmap=cm.RdYlBu_r)
    cbar =plt.colorbar()
    cbar.set_label('({:s})'.format(vari.units.strip()))
    #plt.contour(xx/1e3, yy/1e3, PW_tave, levels=[moist_edgePW], colors='k')
    #plt.plot(mcenter[0], mcenter[1], 'x', markersize=20, zorder=2)
    #plt.title('{:s} ({:s}), day {:3.0f} to {:3.0f} average'.format(varname, vari.units.strip(), t2D[0], t2D[-1]))
    plt.title('{:s} ({:s}), day {:3.0f} to {:3.0f} average'.format(varname, vari.units.strip(), t2D[0], t2D[-1]))
    plt.xlabel('x (km)')
    plt.ylabel('y (km)')
    plt.savefig(fout + '{:s}_day{:3.0f}to{:3.0f}.pdf'.format(varname, t2D[0], t2D[-1]))
    plt.close()

    plt.figure(2)
    ax=plt.gcf().gca()
    if varname == 'PW':
        ax.set_ylim([0,80])
    if varname == 'USFC':
        ax.set_ylim([0,6])
    if varname == 'W500':
        ax.set_ylim([-0.05, 0.3])
        #ax.set_ylim([-0.006, 0.006])
    if varname == 'LHF':
        ax.set_ylim([60,180])
    if varname == 'ZC':
        ax.set_ylim([0, 14])
    if varname == 'Prec':
        ax.set_ylim([0, 300])
    if varname == 'CWP':
        ax.set_ylim([0, 1.2])
    #plt.plot(fieldrs/1e3, fieldprof, 'k,')
    plt.plot(fieldrbin_centers/(domsize*1e3), fieldmeans, 'k.')
    plt.xlabel('fractional distance from moist region center, relative to domain size')
    #plt.xlabel('radial distance (km)')
    plt.ylabel('{:s} ({:s})'.format(varname, vari.units.strip()))
    #plt.title('{:s} ({:s}), day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$, bin width = {:2.2f}'.format(varname, vari.units.strip(), t2D[0], t2D[-1], domsize, 1./binwidth))
    plt.title('{:s} ({:s}), day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$, bin width = {:2.2f} km'.format(varname, vari.units.strip(), t2D[0], t2D[-1], domsize, binwidth/1e3))
    if varname == 'W500':
      ##ax.set_yscale('log')
      fieldrbin_centers = fieldrbin_centers[~np.isnan(fieldmeans)]
      fieldmeans = fieldmeans[~np.isnan(fieldmeans)]
      popt, pcov = curve_fit(expfunc, fieldrbin_centers/(domsize*1e3), fieldmeans)
      a=popt[0]
      b=popt[1]
      c=popt[2]
      fieldfit = expfunc(fieldrbin_centers/(domsize*1e3), a, b, c)
      #fieldfit = expfunc(fieldrbin_centers, a, b, c)
      #plt.plot(fieldrbin_centers/(domsize*1e3), fieldfit, 'b-', alpha=0.6, label=r'${:2.3f}e^{{-{:2.3f}x}} + ({:2.4f})$'.format(a,b,c))
      plt.plot(fieldrbin_centers/(domsize*1e3), fieldfit, 'b-', alpha=0.6, label=r'${:2.3f}e^{{-{:2.3f}x}} + ({:2.4f})$'.format(a,b,c))
      plt.legend()
    plt.savefig(fout + '{:s}radialprof_day{:3.0f}to{:3.0f}.pdf'.format(varname, t2D[0], t2D[-1]))
    ax.set_xlim([0, (1/np.sqrt(2))])
    plt.savefig(fout + '{:s}radialprof_day{:3.0f}to{:3.0f}.pdf'.format(varname, t2D[0], t2D[-1]))
    plt.close()

#plot PW    
#fig = plt.figure()
#ax = fig.gca()
#plt.contourf(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, PW_blocked, 60, cmap=cm.RdYlBu_r)
#cbar = plt.colorbar()
#cbar.set_label('(mm)') 
#plt.contour(xx/1e3, yy/1e3, PW_tave, levels=[moist_edgePW], colors='k')
##plt.plot(mcenter[0], mcenter[1], 'x', markersize=20, zorder=2)
#plt.xlim([x[0], x[-1]])
#plt.ylim([y[0], y[-1]])
##plt.title('PW (mm), day {:3.0f} to {:3.0f} average'.format(t2D[0], t2D[-1]))
#plt.title('PW (mm), day {:3.0f} to {:3.0f} average'.format(t2D[-1], t2D[-1]))
#plt.xlabel('x (km)')
#plt.ylabel('y (km)')
#plt.savefig(fout + 'PW_day{:3.0f}to{:3.0f}.pdf'.format(t2D[0], t2D[-1]))
#plt.close()

#plt.figure(3)
#plt.contourf(xx/1e3, yy/1e3, d/1e3, 60, cmap=cm.RdYlBu_r)
#plt.title('distances from moist region center')
#plt.show()

#plt.figure()
##plt.plot(PWrs/1e3, PWrprof, 'k,')
#plt.plot(PWrbin_centers/(domsize*1e3), PWmeans, 'k.')
#plt.xlabel('fractional distance from moist region center, relative to domain size')
#plt.ylabel('PW (mm)')
#plt.title('PW (mm), day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$, bin width = {:2.2f}'.format(t2D[0], t2D[-1], domsize, 1./binwidth))
#plt.savefig(fout + 'PWradialprof_day{:3.0f}to{:3.0f}.pdf'.format(t2D[0], t2D[-1]))
#plt.close()
#
##plot QNTOA. net radiation at TOA
#fig = plt.figure()
#ax = fig.gca()
#plt.contourf(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, QNTOA_blocked, 60, cmap=cm.RdYlBu_r)
#cbar = plt.colorbar()
#cbar.set_label('(W/m2)')
#plt.contour(xx/1e3, yy/1e3, PW_tave, levels=[moist_edgePW], colors='k')
##plt.plot(mcenter[0], mcenter[1], 'x', markersize=20, zorder=2)
#plt.title('QNTOA (W/m2), day {:3.0f} to {:3.0f} average'.format(t2D[0], t2D[-1]))
#plt.xlabel('x (km)')
#plt.ylabel('y (km)')
#plt.savefig(fout + 'QNTOA_day{:3.0f}to{:3.0f}.pdf'.format(t2D[0], t2D[-1]))
#plt.close()

plt.figure()
plt.plot(QNTOArbin_centers/(domsize*1e3), QNTOAmeans, 'k.')
plt.xlabel('fractional distance from moist region center, relative to domain size')
plt.ylabel('QNTOA (W/m2)')
plt.title('QNTOA (W/m2), day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$, bin width = {:2.2f} km'.format(t2D[0], t2D[-1], domsize, binwidth/1e3))
plt.savefig(fout + 'QNTOAradialprof_day{:3.0f}to{:3.0f}.pdf'.format(t2D[0], t2D[-1]))
plt.close()

fig = plt.figure()
ax = fig.gca()
plt.contourf(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, QN_blocked, 60, cmap=cm.RdYlBu_r)
cbar = plt.colorbar()
cbar.set_label('(W/m2)')
#plt.contour(xx/1e3, yy/1e3, PW_tave, levels=[moist_edgePW], colors='k')
#plt.plot(mcenter[0], mcenter[1], 'x', markersize=20, zorder=2)
plt.title('QNet (W/m2), day {:3.0f} to {:3.0f} average'.format(t2D[0], t2D[-1]))
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.savefig(fout + 'QNet_day{:3.0f}to{:3.0f}.pdf'.format(t2D[0], t2D[-1]))
plt.close()

plt.figure()
plt.plot(QNrbin_centers/(domsize*1e3), QNmeans, 'k.')
plt.xlabel('fractional distance from moist region center, relative to domain size')
plt.ylabel('QNet (W/m2)')
plt.title('QNet (W/m2), day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$, bin width = {:2.2f} km'.format(t2D[0], t2D[-1], domsize, binwidth/1e3))
plt.savefig(fout + 'QNetradialprof_day{:3.0f}to{:3.0f}.pdf'.format(t2D[0], t2D[-1]))
plt.close()

fig = plt.figure()
ax = fig.gca()
plt.contourf(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, SWN_blocked, 60, cmap=cm.RdYlBu_r)
cbar = plt.colorbar()
cbar.set_label('(W/m2)')
plt.contour(xx/1e3, yy/1e3, PW_tave, levels=[moist_edgePW], colors='k')
#plt.plot(mcenter[0], mcenter[1], 'x', markersize=20, zorder=2)
plt.title('SWN (W/m2), day {:3.0f} to {:3.0f} average'.format(t2D[0], t2D[-1]))
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.savefig(fout + 'SWN_day{:3.0f}to{:3.0f}.pdf'.format(t2D[0], t2D[-1]))
plt.close()

plt.figure()
plt.plot(QNrbin_centers/(domsize*1e3), SWNmeans, 'k.')
plt.xlabel('fractional distance from moist region center, relative to domain size')
plt.ylabel('SWN (W/m2)')
plt.title('SWN (W/m2), day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$, bin width = {:2.2f} km'.format(t2D[0], t2D[-1], domsize, binwidth/1e3))
plt.savefig(fout + 'SWNradialprof_day{:3.0f}to{:3.0f}.pdf'.format(t2D[0], t2D[-1]))
plt.close()








