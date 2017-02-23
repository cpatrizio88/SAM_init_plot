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

matplotlib.rcParams.update({'font.size': 24})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 2})
matplotlib.rcParams.update({'legend.fontsize': 22})

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fpath3D = '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpathSTAT = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
foutdata = '/Users/cpatrizio/data/SST302/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_RADIAL/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to140_64vert_ubarzero_RADIAL/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to150_64vert_ubarzero_RADIAL/'

fout='/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_RADIAL/'



#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

nc_inSTAT1 = glob.glob(fpathSTAT + '*256x256*3000m*250days*302K.nc')[0]
nc_in2D1 = glob.glob(fpath + '*256x256*3000m*day230to250*302K.nc')[0]
nc_in3D1 = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]

nc_inSTAT2 = glob.glob(fpathSTAT + '*512x512*3000m*180days*302K.nc')[0]
nc_in2D2 = glob.glob(fpath + '*512x512*3000m*day180to195*302K.nc')[0]
nc_in3D2 = glob.glob(fpath3D + '*512x512*3000m*day180to195*302K.nc')[0]

nc_inSTAT3 = glob.glob(fpathSTAT + '*1024x1024*3000m*180days*302K.nc')[0]
nc_in2D3 = glob.glob(fpath + '*1024x1024*3000m*day190to200*302K.nc')[0]
nc_in3D3 = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]

domsizes = [768, 1536, 3072]
nc_fs = [nc_in2D1, nc_in2D2, nc_in2D3]
nc_fs3D = [nc_in3D1, nc_in3D2, nc_in3D3]
nc_STATs = [nc_inSTAT1, nc_inSTAT2, nc_inSTAT3]
colors = ['k', 'r', 'g']

varnames = ['PW', 'USFC', 'LHF', 'CWP','W500', 'TB', 'ZC', 'Prec']
#varnames = ['ZC']
#varnames = ['UBLRAD']
#varnames = ['PSFC']
#varnames = ['USFC']

#varnames = ['W']

#varnames = ['LWNTC', 'LWNSC', 'SWNTC', 'SWNSC']

#varnames = ['LWNS', 'LWNT', 'SWNS', 'SWNT']

#varnames = ['Prec']

colors = ['k', 'r', 'g']

for i, domsize in enumerate(domsizes): 
    
    print 'domsize', domsize
    
    print 'loading .nc file', nc_fs[i]
    print 'loading 2D variables'
    nc_data2D = Dataset(nc_fs[i])
    varis2D = nc_data2D.variables
    
    print 'loading 3D variables'
    nc_data3D = Dataset(nc_fs3D[i])
    varis3D = nc_data3D.variables
    
    nc_inSTAT = nc_STATs[i]
    nc_dataSTAT = Dataset(nc_inSTAT)
    varisSTAT = nc_dataSTAT.variables
    
    z = varisSTAT['z'][:]
    p = varisSTAT['p'][:]
    p = p*1e2
    
    if domsize == 768:
        nave=10
    elif domsize == 1536:
        nave=10
    else:
        nave=3
        
    ntave2D=24
    ntave3D=4
    t=-1
    t3=-1

    
    aveperiod = nave*ntave2D
    aveperiod3D = nave*ntave3D
    
    PW = varis2D['PW'][t-aveperiod:t,:,:]
    USFC = varis2D['USFC'][t-aveperiod:,:,:]
    VSFC = varis2D['VSFC'][t-aveperiod:,:,:]
    
    speedSFC= np.sqrt(np.power(USFC, 2) + np.power(VSFC,2))
    LWNT = varis2D['LWNT'][t-aveperiod:t,:,:]
    LWNS = varis2D['LWNS'][t-aveperiod:t,:,:]
    SWNT = varis2D['SWNT'][t-aveperiod:t,:,:]
    SWNS = varis2D['SWNS'][t-aveperiod:t,:,:]
    
    LWNTC = varis2D['LWNTC'][t-aveperiod:t,:,:]
    LWNSC = varis2D['LWNSC'][t-aveperiod:t,:,:]
    
    SWNTC = varis2D['SWNTC'][t-aveperiod:t,:,:]
    SWNSC = varis2D['SWNSC'][t-aveperiod:t,:,:]
        
    QNTOA = SWNT - LWNT
    QNS = SWNS - LWNS
    
    SWN = SWNT - SWNS
    LWN = LWNS - LWNT
    
    LWNC = LWNSC - LWNTC
    SWNC = SWNTC - SWNSC
    
    QN = QNTOA - QNS
        
    #field = np.sqrt(np.power(USFC, 2) + np.power(VSFC, 2))
    t2D = varis2D['time'][t-aveperiod:t]
    x = varis2D['x'][:]
    y = varis2D['y'][:]
    nt2D = t2D.size
    nx = x.size
    ny = y.size
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
    QNTOA_tave = np.mean(QNTOA, axis=0)
    QN_tave = np.mean(QN, axis=0)
    SWN_tave = np.mean(SWN, axis=0)
    LWN_tave = np.mean(LWN, axis=0)
    SWNC_tave = np.mean(SWNC, axis=0)
    LWNC_tave = np.mean(LWNC, axis=0)
    
    xx, yy = np.meshgrid(x, y)
    times = np.arange(t2D[0], t2D[-1])
    db=1
    #AVERAGING PERIOD IN DAYS
    
    #PW_t = PW_tave[t,:,:]
    #PW_t = np.mean(PW_tave[t-nave:t,:,:], axis=0)
    #QNTOA_t = QNTOA_tave[t,:,:]
    #QNTOA_t = np.mean(QNTOA_tave[t-nave:t,:,:], axis=0)
    #QN_t = np.mean(QN_tave[t-nave:t,:,:], axis=0)
    
    PW_blocked = blockave2D(PW_tave, db)
    QNTOA_blocked=blockave2D(QNTOA_tave, db)
    SWN_blocked = blockave2D(SWN_tave, db)
    LWN_blocked = blockave2D(LWN_tave, db)
    QN_blocked=blockave2D(QN_tave, db)
    SWNC_blocked = blockave2D(SWNC_tave, db)
    LWNC_blocked = blockave2D(LWNC_tave, db)
    
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
    PWrbins, PWmeans = radprof(PW_tave, xx, yy, mcenter, bins)
    PWrbin_centers = (PWrbins[1:] + PWrbins[:-1])/2.
    #PWrprof = PWrsums/PWnr
    #PWrs = np.arange(0, len(PWnr))
    
    QNTOArbins, QNTOAmeans = radprof(QNTOA_tave, xx, yy, mcenter, bins)
    QNTOArbin_centers = (QNTOArbins[1:] + QNTOArbins[:-1])/2.
    #QNTOArprof = QNTOAsums/QNTOAnr
    #QNTOArs = np.arange(0, len(QNTOAnr))
    
    QNrbins, QNmeans = radprof(QN_tave, xx, yy, mcenter, bins)
    QNrbin_centers = (QNrbins[1:] + QNrbins[:-1])/2.
    
    SWNrbins, SWNmeans = radprof(SWN_tave, xx, yy, mcenter, bins)
    SWNrbin_centers = (SWNrbins[1:] + SWNrbins[:-1])/2.
    
    LWNrbins, LWNmeans = radprof(LWN_tave, xx, yy, mcenter, bins)
    LWNrbin_centers = (LWNrbins[1:] + LWNrbins[:-1])/2.
    
    LWNCrbins, LWNCmeans = radprof(LWNC_tave, xx, yy, mcenter, bins)
    LWNCrbin_centers = (LWNCrbins[1:] + LWNCrbins[:-1])/2.
    
    SWNCrbins, SWNCmeans = radprof(SWNC_tave, xx, yy, mcenter, bins)
    SWNCrbin_centers = (SWNCrbins[1:] + SWNCrbins[:-1])/2.
    
    a=1.5 
    moist_edgePW = np.mean(PW_tave) + a*np.std(PW_tave)  
    
    print 'calculating radial profiles'
    
    for vi, varname in enumerate(varnames):

        print varname
        if varname == 'W':
            
            
            vari = varis3D[varname]
            units = vari.units.strip()
            field = vari[t3-aveperiod3D:t3,:,:,:]
            field_tave = np.mean(field, axis=0)
            #take vertical average
            #field_ztave = np.mean(field_tave, axis=0)
            
            #take free troposphere average?
            p_t = 150
            z_BL = 700
            BLi = np.where(z > z_BL)[0][0]
            #BLi=0
            ti = np.where(p <= p_t*1e2)[0][0]
            #ti=-1
            field_tave = np.mean(field_tave[BLi:ti,:,:], axis=0)
            #field_blocked = blockave3D(field_tave, db)

        #calculate speed instead of zonal wind velocity at surface
        elif varname == 'USFC': 
            field=speedSFC
            vari = varis2D['USFC']
            field_tave = np.mean(field, axis=0)
        elif varname == 'UBLRAD':
            U_rfname = glob.glob(foutdata + '*{:d}*_URADIAL_*{:3.0f}to{:3.0f}*'.format(domsize, t2D[0], t2D[-1]))[0]
            vari = varis2D['USFC']
            print 'loading', U_rfname
            field_tave = np.load(U_rfname)   
            z_BL = 500        
            BLi = np.where(z > z_BL)[0][0]
            #field_tave = np.mean(field_tave, axis=0)
            field_tave = np.mean(field_tave[:BLi,:,:], axis=0)
            #field_tave = blockave3D(field_tave, db)
            
        else:
            vari = varis2D[varname]
            field = varis2D[varname][t-aveperiod:t,:,:]
            field_tave = np.mean(field, axis=0)
        
        #field_t = field_tave[t,:,:]
        #average over nave days 
        #field_t = np.mean(field_tave[t-nave:t,:,:], axis=0)
        
        #field_blocked = blockave2D(field_tave, db)
        
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
        #elif varname == 'W':
        #    cvals = np.arange(-0.003, 0, 0.0001)
        else:
            cvals = 60
        
        
        #fig = plt.figure(1)
        #ax = fig.gca()
        #plt.contourf(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, field_blocked, cvals, cmap=cm.RdYlBu_r)
        #cbar =plt.colorbar()
        #cbar.set_label('({:s})'.format(vari.units.strip()))
        ##plt.contour(xx/1e3, yy/1e3, PW_tave, levels=[moist_edgePW], colors='k')
        ##plt.plot(mcenter[0], mcenter[1], 'x', markersize=20, zorder=2)
        ##plt.title('{:s} ({:s}), day {:3.0f} to {:3.0f} average'.format(varname, vari.units.strip(), t2D[0], t2D[-1]))
        #plt.title('{:s} ({:s}), day {:3.0f} to {:3.0f} average'.format(varname, vari.units.strip(), t2D[0], t2D[-1]))
        #plt.xlabel('x (km)')
        #plt.ylabel('y (km)')
        #plt.savefig(fout + '{:s}_day{:3.0f}to{:3.0f}.pdf'.format(varname, t2D[0], t2D[-1]))
        #plt.close()
    
        plt.figure(vi)
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
        if varname == 'W':
            ax.set_ylim([-0.05, 0.3])
            varname = 'tropospheric average W'
        if varname == 'UBLRAD':
            titlename = r'u$_{BL,r}$'
        else: 
            titlename = varname
            #ax.set_ylim([-0.006, 0])
        #plt.plot(fieldrs/1e3, fieldprof, 'k,')
        plt.plot(fieldrbin_centers/(domsize*1e3), fieldmeans, 'x', mew=2, color=colors[i],  label='{:d} km, day {:3.0f} to {:3.0f} mean'.format(domsize, t2D[0], t2D[-1]))
        plt.xlabel(r'$\hat{r}$')
        #plt.xlabel('radial distance (km)')
        plt.ylabel('{:s} ({:s})'.format(titlename, vari.units.strip()))
        #plt.title('{:s} ({:s}), day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$, bin width = {:2.2f}'.format(varname, vari.units.strip(), t2D[0], t2D[-1], domsize, 1./binwidth))
        plt.title('{:s}'.format(titlename))
        #if varname == 'W500':
        ##ax.set_yscale('log')
        #fieldrbin_centers = fieldrbin_centers[~np.isnan(fieldmeans)]
        #fieldmeans = fieldmeans[~np.isnan(fieldmeans)]
        #popt, pcov = curve_fit(expfunc, fieldrbin_centers/(domsize*1e3), fieldmeans)
        #a=popt[0]
        #b=popt[1]
        #c=popt[2]
        #fieldfit = expfunc(fieldrbin_centers/(domsize*1e3), a, b, c)
        #fieldfit = expfunc(fieldrbin_centers, a, b, c)
        #plt.plot(fieldrbin_centers/(domsize*1e3), fieldfit, 'b-', alpha=0.6, label=r'${:2.3f}e^{{-{:2.3f}x}} + ({:2.4f})$'.format(a,b,c))
        #plt.plot(fieldrbin_centers/(domsize*1e3), fieldfit, 'b-', alpha=0.6, label=r'${:2.3f}e^{{-{:2.3f}x}} + ({:2.4f})$'.format(a,b,c))
        #plt.legend()
        #plt.savefig(fout + '{:s}radialprof_day{:3.0f}to{:3.0f}.pdf'.format(varname, t2D[0], t2D[-1]))
        ax.set_xlim([0, (1/np.sqrt(2))])
        
        plt.savefig(fout + '{:s}radialprof_day250_{:d}day.pdf'.format(varname, nave))
        #plt.close()
    
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
    
    plt.figure(16)
    plt.plot(QNTOArbin_centers/(domsize*1e3), QNTOAmeans, '-', color=colors[i],  label='{:d} km, day {:3.0f} to {:3.0f} mean'.format(domsize, t2D[0], t2D[-1]))
    plt.xlabel('fractional distance from moist region center, relative to domain size')
    plt.ylabel('QNTOA (W/m2)')
    plt.title('QNTOA (W/m2), bin width = {:2.2f} km'.format(binwidth/1e3))
    plt.savefig(fout + 'QNTOAradialprof_day250.pdf_{:d}day.pdf'.format(nave))

    
    #fig = plt.figure(11)
    #ax = fig.gca()
    #plt.contourf(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, QN_blocked, 60, cmap=cm.RdYlBu_r)
    #cbar = plt.colorbar()
    #cbar.set_label('(W/m2)')
    ##plt.contour(xx/1e3, yy/1e3, PW_tave, levels=[moist_edgePW], colors='k')
    ##plt.plot(mcenter[0], mcenter[1], 'x', markersize=20, zorder=2)
    #plt.title('QNet (W/m2), day {:3.0f} to {:3.0f} average'.format(t2D[0], t2D[-1]))
    #plt.xlabel('x (km)')
    #plt.ylabel('y (km)')
    #plt.savefig(fout + 'QNet_day{:3.0f}to{:3.0f}.pdf'.format(t2D[0], t2D[-1]))
    #plt.close()
    
    plt.figure(17)
    plt.plot(QNrbin_centers/(domsize*1e3), QNmeans, '-', color=colors[i],  label='{:d} km, day {:3.0f} to {:3.0f} mean'.format(domsize, t2D[0], t2D[-1]))
    plt.xlabel(r'$\hat{r}$')
    plt.ylabel(r'$Q_{r,net}$ (W/m$^2$)')
    plt.title(r'net radiative heating')
    plt.savefig(fout + 'QNetradialprof_day250_{:d}day.pdf'.format(nave))

    
    #fig = plt.figure()
    #ax = fig.gca()
    #plt.contourf(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, SWN_blocked, 60, cmap=cm.RdYlBu_r)
    #cbar = plt.colorbar()
    #cbar.set_label('(W/m2)')
    #plt.contour(xx/1e3, yy/1e3, PW_tave, levels=[moist_edgePW], colors='k')
    ##plt.plot(mcenter[0], mcenter[1], 'x', markersize=20, zorder=2)
    #plt.title('SWN (W/m2), day {:3.0f} to {:3.0f} average'.format(t2D[0], t2D[-1]))
    #plt.xlabel('x (km)')
    #plt.ylabel('y (km)')
    #plt.savefig(fout + 'SWN_day{:3.0f}to{:3.0f}.pdf'.format(t2D[0], t2D[-1]))
    #plt.close()
    
    plt.figure(18)
    plt.plot(QNrbin_centers/(domsize*1e3), SWNmeans, '-', color=colors[i],  label='{:d} km, day {:3.0f} to {:3.0f} mean'.format(domsize, t2D[0], t2D[-1]))
    plt.xlabel(r'$\hat{r}$')
    plt.ylabel(r'$Q_{SW,net}$ (W/m$^2$)')
    plt.title('net shortwave heating')
    plt.savefig(fout + 'SWNradialprof_day250.pdf_{:d}day.pdf'.format(nave))
    
    
    plt.figure(19)
    plt.plot(QNrbin_centers/(domsize*1e3), LWNmeans, '-', color=colors[i],  label='{:d} km, day {:3.0f} to {:3.0f} mean'.format(domsize, t2D[0], t2D[-1]))
    plt.xlabel(r'$\hat{r}$')
    plt.ylabel(r'$Q_{LW,net}$ (W/m$^2$)')
    plt.title('net longwave heating')
    plt.savefig(fout + 'LWNradialprof_day250.pdf_{:d}day.pdf'.format(nave))
    
    plt.figure(20)
    plt.plot(QNrbin_centers/(domsize*1e3), SWNCmeans, '-', color=colors[i],  label='{:d} km, day {:3.0f} to {:3.0f} mean'.format(domsize, t2D[0], t2D[-1]))
    plt.xlabel(r'$\hat{r}$')
    plt.ylabel(r'$Q_{SWcl,net}$ (W/m$^2$)')
    plt.title('cloud net shortwave heating')
    plt.savefig(fout + 'SWNCradialprof_day250.pdf_{:d}day.pdf'.format(nave))

    
    plt.figure(21)
    plt.plot(QNrbin_centers/(domsize*1e3), LWNCmeans, '-', color=colors[i],  label='{:d} km, day {:3.0f} to {:3.0f} mean'.format(domsize, t2D[0], t2D[-1]))
    plt.xlabel(r'$\hat{r}$')
    plt.ylabel(r'$Q_{LWcl,net}$ (W/m$^2$)')
    plt.title('cloud net longwave heating')
    plt.savefig(fout + 'LWNCradialprof_day250.pdf_{:d}day.pdf'.format(nave))
    
    
    plt.figure(22)
    plt.plot(QNrbin_centers/(domsize*1e3), LWNmeans-LWNCmeans, '-', color=colors[i],  label='{:d} km, day {:3.0f} to {:3.0f} mean'.format(domsize, t2D[0], t2D[-1]))
    plt.xlabel(r'$\hat{r}$')
    plt.ylabel(r'$Q_{LWclear,net}$ (W/m$^2$)')
    plt.title('clear-sky net longwave heating')
    plt.savefig(fout + 'LWNclearradialprof_day250_{:d}day.pdf'.format(nave))
    
    
    plt.figure(23)
    plt.plot(QNrbin_centers/(domsize*1e3), SWNmeans-SWNCmeans, '-', color=colors[i],  label='{:d} km, day {:3.0f} to {:3.0f} mean'.format(domsize, t2D[0], t2D[-1]))
    plt.xlabel('fractional distance from moist region center, relative to domain size')
    plt.ylabel('clear-sky net shortwave heating (W/m2)')
    plt.title('clear-sky shortwave heating (W/m2), bin width = {:2.2f} km'.format(binwidth/1e3))
    plt.savefig(fout + 'SWNclearradialprof_day250_{:d}day.pdf'.format(nave))
    
    
    
for vi, varname in enumerate(varnames):
    plt.figure(vi)
    plt.legend(loc='best')
    plt.savefig(fout + '{:s}radialprof_day250_{:d}day.pdf'.format(varname, nave))
    plt.close()
    
varnames = ['QNTOA', 'QNet', 'SWN', 'LWN', 'SWNC', 'LWNC', 'LWNclear', 'SWNclear']
    
for i, fi in enumerate(np.arange(16,16+len(varnames))):
    plt.figure(fi)
    plt.legend(loc='best')
    plt.savefig(fout + '{:s}radialprof_day250_{:d}day.pdf'.format(varnames[i], nave))
    plt.close()








