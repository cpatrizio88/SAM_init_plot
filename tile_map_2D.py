from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker
import matplotlib
from SAM_init_plot.block_fns import blockave3D, blockave2D
from thermolib.constants import constants

from matplotlib.ticker import AutoMinorLocator

c = constants()

matplotlib.rcParams.update({'font.size': 26})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 2})
matplotlib.rcParams.update({'legend.fontsize': 22})

plt.style.use('seaborn-white')
fpathSTAT = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
fpath2D =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fpath3D = '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
foutdata = '/Users/cpatrizio/data/SST302/'
#fout =  '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_MAPSNEW/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to140_64vert_ubarzero_MAPSNEW/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to150_64vert_ubarzero_MAPSNEW/'


#nc_STAT = glob.glob(fpathSTAT + '*256x256*3000m*250days*302K.nc')[0]
nc_in1 = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]
 
#nc_STAT = glob.glob(fpathSTAT + '*512x512*3000m*195days*302K.nc')[0]
nc_in2 = glob.glob(fpath2D + '*512x512*3000m*day180to195*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day090to130*302K.nc')[0]

#nc_STAT = glob.glob(fpathSTAT + '*1024x1024*3000m*220days*302K.nc')[0]
nc_in3 = glob.glob(fpath2D + '*1024x1024*3000m*day170to180*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]

#domsize=768
#domsize=1536
#domsize=3072

nc_data1= Dataset(nc_in1)
nc_data2 = Dataset(nc_in2)
nc_data3 = Dataset(nc_in3)
#nc_dataSTAT = Dataset(nc_STAT)
#nc_data3D = Dataset(nc_in3D)
varis2D1 = nc_data1.variables
varis2D2 = nc_data2.variables
varis2D3 = nc_data3.variables
#varis3D = nc_data3D.variables
#varisSTAT = nc_dataSTAT.variables
#times2D = varis2D['time'][:]

#times3D = varis3D['time'][:]

times2D2 = varis2D2['time'][:]

x1 = varis2D1['x'][:]
y1 = varis2D1['y'][:]

x2 = varis2D2['x'][:]
y2 = varis2D2['y'][:]

x3 = varis2D3['x'][:]
y3 = varis2D3['y'][:]


#fac=12

#nave2D = 24/fac
#nave3D = 4/fac
nave2D=24
nave3D=4

#times = np.arange(nave2D, 15*fac)

db=1

t1=90*nave2D
t2=-1
t3=-1


xx, yy = np.meshgrid(x3, y3)
#LWNT: OLR
#Prec: surface precip
#varname = 'W500'
#vari = varis[varname]
#field = vari[:]
#field_tave = blockave3D(field, db)

#
#p = varis3D['p'][:]
#z = varis3D['z'][:]
varname = 'SFCSPEED'
varname = 'PW'
vari = varis2D1[varname]
times2D1 = varis2D1['time'][:]
times2D3 = varis2D3['time'][:]
#plev = 500
#plevi = np.where(p > plev)[0][-1]
#wfield_p = varis3D['W'][:,plevi,:,:]
#block_field = blockave3D(wfield_p, db)
#field = block_field

#varname = 'NETSW'


nx=x3.size
ny=y3.size
#nz=z.size


#tstep=int(4)
#tstep=int(12)
#tstart=4*90
#tstart = 24*40
#tstart=nave3D


#W_crit = 0.01
#Wedge_u = 3
#Wedge_l = 0.5
#vals=[Wedge_l*W_crit, 0.01, Wedge_u*W_crit]
#vals=20

#varname = 'SFCSPEED'
#varname = 'BLQVADV'

#CAPE COMPUTATION
#varname = 'CAPE'
#TABS = varis3D['TABS'][:,:,:,:]
#QV = varis3D['QV'][:,:,:,:]
#QV = QV*1e-3
#
#TV = TABS*(1+0.61*QV)

#varname ='URADIAL'
#varname = 'Nsquared'
if varname == 'CAPE':
      units = r'J/m$^2$'
elif varname == 'N':
      units = r's$^{{-2}}$'
elif varname == 'URADIAL':
      U_rfname = glob.glob(foutdata + '{:d}km_URADIAL_*{:3.0f}*'.format(domsize, times3D[-1]))[0] 
      U_r = np.load(U_rfname) 
      units = 'm/s'
elif varname == 'SFCSPEED':
      units = 'm/s'
elif varname == 'BLQVADV':
      units = r'g kg$^{-1}$ s$^{-1}$'
elif varname == 'NETLW':
      units = r'W/m$^{2}$'
elif varname == 'NETSW':
      units = r'W/m$^{2}$'
elif varname == 'ALPHA':
      units=''
else:
      units = vari.units.strip()

#delz = np.diff(z)

#vals = np.linspace(-2, 5, 50)

#PW values
vals = np.arange(0, 88, 4)

#CWP values
#vals = np.arange(0, 1.3, 0.02)

#N^2 values
#vals = np.linspace(2.3e-4, 2.6e-4, 50)

#vals = np.linspace(-1, 1, 50)

#ZC vals
#vals = np.arange(0, 15, 0.2)


#CAPE vals?
#vals = np.arange(0, 800, 10)

#Prec vals
#vals = np.arange(0, 1200, 20)

#W blockave vals
#vals = np.arange(-0.10, 0.5, 0.01)

#W vals
#vals = np.arange(-0.5, 8, 0.1)

#USFC vals
#vals= np.arange(0, 9, 0.2)

#ALPHA VALS
#vals = np.linspace(0, 0.1, 50)


#NETLW vals
#vals=np.arange(-290,-50, 10)

#vals = 50

#vals = np.arange(-0.001, 0.001, .0001)

#alpha_t = np.zeros(times.size)

#for ii, i in enumerate(np.arange(tstart, len(times3D), tstep)):
#EDIT THIS TO PLOT PROJECTION OF 3D FIELD
#for i in times:
#    
#    t2=i*nave2D
#    print times2D[t2]
#    t3=i*nave3D
#    
#    if varname == 'CAPE':
#     
#        #CAPE CALCULATION
#        TABS = varis3D['TABS'][i-nave3D:i,:,:,:]
#        QV = varis3D['QV'][i-nave3D:i,:,:,:]
#        QV = QV*1e-3
#        
#        TABS_tave = np.mean(TABS, axis=0)
#        QV_tave = np.mean(QV, axis=0)
#    
#        TV_tave = TABS_tave*(1+0.61*QV_tave)
#    
#        TENV = nc_STAT['TABS'][i-nave2D:i,:]
#        TENV = np.mean(TENV, axis=0)
#        QBAR = nc_STAT['QV'][i-nave2D:i,:]
#        QBAR = np.mean(QBAR, axis=0)
#        TVENV = TENV*(1+0.61*QBAR)
#        TVENV3D = np.zeros((ny, nx, nz))
#        TVENV3D[:,:,:] = TVENV
#        TVENV3D = TVENV3D.T
#        
#        CAPE = c.g*np.sum(np.multiply((TV_tave - TVENV)/TVENV, delz))
#        
#        field_tave = CAPE
#     
#    elif varname == 'Nsquared':
#        p0=1000
#        T = varis3D['TABS'][t3-nave3D:t3,:,:,:]
#        T_tave = np.mean(T, axis=0)
#        p3D = np.zeros((nz, nx, ny))
#        p3D= p3D.T
#        p3D[:,:,:] = p
#        p3D = p3D.T
#        #p3D = p3D.T
#        #p4D = np.tile(p3D[:,:,:, np.newaxis], T.shape[0])
#        #p4D = p4D.T
#        THETA_tave = T_tave*(p0/p3D)**(c.Rd/c.cp)
#        diffTHETA = THETA_tave[1:,:,:] - THETA_tave[:-1,:,:]
#        diffz = np.diff(z)
#        diffz3D = np.zeros((nz-1, nx, ny))
#        diffz3D = diffz3D.T
#        diffz3D[:,:,:] = diffz
#        diffz3D = diffz3D.T
#        THETA_mid = (THETA_tave[1:,:,:] + THETA_tave[:-1,:,:])/2.
#        N_tave = (c.g*diffTHETA)/np.multiply(diffz3D, THETA_mid)
#        N_tave_vertave = np.mean(N_tave, axis=0)
#        field_tave = blockave2D(N_tave_vertave, db)
#        
#    elif varname == 'URADIAL':
#        U_r_p = U_r[plevi,:,:]
#        field_tave = blockave2D(U_r_p, db)
#        
#    
#    elif varname == 'SFCSPEED':
#        USFC = varis2D['USFC'][t2-nave2D:t2,:,:]
#        VSFC = varis2D['VSFC'][t2-nave2D:t2,:,:]
#        speed = np.sqrt(np.power(USFC, 2) + np.power(VSFC, 2))
#        field_tave = np.mean(speed, axis=0)
#        field_tave = blockave2D(field_tave, db)
#        
#        #look at different levels here 
#        #p = varis3D['p'][:]
#        #z = varis3D['z'][:]
#        #varname = 'W'#
#        #plev = 500
#        #plevi = np.where(p > plev)[0][-1]
#        #U = varis3D['U'][t3-nave3D:t3,:,:,:]
#        #V = varis3D['V'][t3-nave3D:t3,:,:,:]
#        
#    elif varname == 'BLQVADV':
#        QV = varis3D['QV'][t3-nave3D:t3,:,:,:]
#        U = varis3D['U'][t3-nave3D:t3,:,:,:]
#        V = varis3D['V'][t3-nave3D:t3,:,:,:]
#        z_BL = 1e3
#        BLi = np.where(z > z_BL)[0][0]
#        QV_BL = np.mean(QV[:,:BLi,:,:], axis=1)
#        U_BL = varis3D['U'][t3-nave3D:t3,:BLi,:,:]
#        V_BL = varis3D['V'][t3-nave3D:t3,:BLi,:,:]
#        U_BL = np.mean(U_BL, axis=1)
#        V_BL = np.mean(V_BL, axis=1)
#        
#        delx = np.diff(x)[0]
#        dely = np.diff(y)[0]
#        
#        #calculate horizontal advection of moisture use forward difference
#        
#        #diffU = U[:,:,1:,:] - U[:,:,:-1,:]
#        #diffU = diffU[:,:,:,:-1]
#        #diffV = V[:,:,:,1:] - V[:,:,:,:-1]
#        #diffV = diffV[:,:,:-1,:]
#        #diffQVx = QV[:,:,1:,:] - QV[:,:,:-1,:]
#        #diffQVx = diffQVx[:,:,:,:-1]
#        #diffQVy = QV[:,:,:,1:] - QV[:,:,:,:-1]
#        #diffQVy =diffQVy[:,:,:-1,:]
#        #QV = QV[:,:,:-1,:-1]
#        
#        diffU = U_BL[:,1:,:] - U_BL[:,:-1,:]
#        diffU = diffU[:,:,:-1]
#        diffV = V_BL[:,:,1:] - V_BL[:,:,:-1]
#        diffV = diffV[:,:-1,:]
#        diffQVx = QV_BL[:,1:,:] - QV_BL[:,:-1,:]
#        diffQVx = diffQVx[:,:,:-1]
#        diffQVy = QV_BL[:,:,1:] - QV_BL[:,:,:-1]
#        diffQVy =diffQVy[:,:-1,:]
#        QV_BL = QV_BL[:,:-1,:-1]
#        
#        
#        U_BL = U_BL[:,:-1,:-1]
#        V_BL = V_BL[:,:-1,:-1]
#        
#        #U = U[:,:,:-1,:-1]
#        #V = V[:,:,:-1,:-1]
#        
#        #QVadv = QV_BL*(diffU/delx + diffV/dely) + U_BL*(diffQVx/delx) + V_BL*(diffQVy/dely)
#        QVadv = U_BL*(diffQVx/delx) + V_BL*(diffQVy/dely)
#        
#        #QVadv = QV*(diffU/delx + diffV/dely) + U*(diffQVx/delx) + V*(diffQVy/dely)
#        #field_tave = np.mean(QVadv, axis=1)
#        #field_tave = np.mean(field_tave, axis=0)
#        field_tave = np.mean(QVadv, axis=0)
#    elif varname == 'NETLW':
#        LWNTOA = varis2D['LWNT'][t2-nave2D:t2,:,:]
#        LWNS = varis2D['LWNS'][t2-nave2D:t2,:,:]
#        netLW = LWNS - LWNTOA
#        field_tave = np.mean(netLW, axis=0) 
#        field_tave = blockave2D(field_tave, db)
#    elif varname == 'NETSW':
#        SWNTOA = varis2D['SWNT'][t2-nave2D:t2,:,:]
#        SWNS = varis2D['SWNS'][t2-nave2D:t2,:,:]
#        netLW = SWNTOA - SWNS
#        field_tave = np.mean(netLW, axis=0)
#        field_tave = blockave2D(field_tave, db)
#    elif varname == 'ALPHA':
#        LHF = varis2D['LHF'][t2-nave2D:t2,:,:]
#        P = varis2D['Prec'][t2-nave2D:t2,:,:]
#        LHF = blockave3D(LHF, db)
#        P = blockave3D(P, db)
#        E = (LHF*1000*86400)/(L_v*rho_w)
#        alpha = E/P
#        field_tave = np.mean(alpha, axis=0)
#        indx = np.isfinite(field_tave)
#        alpha_t[i] = np.mean(field_tave[indx])
#        #field_tave = blockave2D(field_tave, db)
#    else:
#        field = varis2D[varname][t2-nave2D:t2,:,:]
#        field_tave = np.mean(field, axis=0)
#        field_tave = blockave2D(field_tave, db)
#        
#        
#        
#        
#        
#
#        
#        
#
#
#    #w_tave = np.mean(wfield_p[t3-nave3D:t3,:,:], axis=0)
#    #w_tave = blockave2D(w_tave, db)
#    
#    #wvals = [W_crit]

def newtile(field, n):
    old_nx = field.shape[0]
    old_ny = field.shape[1]
    shape = (n*old_nx, n*old_ny)
    new_field = np.zeros(shape)
    for i in range(n):
        for j in range(n):
            new_field[old_nx*j:old_nx*(j+1),old_ny*i:old_ny*(i+1)] = field
            
    return new_field
            
        

PW_field1 = varis2D1['PW'][t1,:,:]
PW_field2 = varis2D2['PW'][t2,:,:]
PW_field3 = varis2D3['PW'][t3,:,:]

PW_field1new = newtile(PW_field1, 4)
#PW_field1new = PW_field1new.reshape(PW_field3.shape)

PW_field2new = newtile(PW_field2, 2)
#PW_field2new = PW_field2new.reshape(PW_field3.shape)

#PW_tave = np.mean(PW_field, axis=0)
#PW_tave = blockave2D(PW_tave, db)

PWvals = np.arange(0, 88, 8)

#major_ticks1 = np.arange(0, 3*1024+256, 256)
major_ticks2 = np.arange(0, 3*1024+768, 768)


#field_tave = np.mean(field[t2-nave2D:t2,:,:], axis=0)
#field_tave = blockave2D(field_tave, db)
#varname_u = 'U200'
#U = varis[varname_u][t2-nave2D:t2,:,:]
#U_tave = np.mean(U, axis=0)
#V = varis['V200'][i-nave2D:i,:,:]
#V_tave = np.mean(V, axis=0)

#print 'day', times2D[t2]
#print 'day', times3D[i]
#plt.figure()
#minvar= np.min(field_tave[~np.isnan(field_tave)])
#maxvar= np.max(field_tave[~np.isnan(field_tave)])
#plt.contour(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, field_tave[:,:], vals, colors='k', linewidth=0.5, alpha=0.5)
titlename='PW'
plt.figure(1)
cs = plt.contourf(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, PW_field1new[:,:], PWvals, cmap=cm.RdYlBu_r, zorder=0)
ax=plt.gca()
#ticklabels = np.arange(0, 3500, 500)
ax.set_xticks(major_ticks2)
ax.set_yticks(major_ticks2)
minorlocator = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(minorlocator)
ax.yaxis.set_minor_locator(minorlocator)
ax.grid(which='major', color='k', zorder=10)
#cs.collections[0].set_label('W$_c$ = {:3.2f} m/s'.format(W_crit))
#cs.collections[0].set_label('PW')
#plt.legend(loc='best')
#q = plt.quiver(xx[::8, ::8]/1e3, yy[::8, ::8]/1e3, U_tave[::8,::8], V_tave[::8,::8], scale=500, alpha=0.8, zorder=1)
if varname == 'BLQVADV':
    xxplot = xx[:-1, :-1]
    yyplot = yy[:-1, :-1]
else:
    xxplot = xx
    yyplot = yy
#plt.contourf(xxplot[::db, ::db]/1e3, yyplot[::db, ::db]/1e3, field_tave[:,:], vals, cmap=cm.RdYlBu_r, zorder=0)
cb = plt.colorbar()
cb.set_label('({0})'.format(units))
#p = plt.quiverkey(q, np.min(x)/1e3+30, np.max(y)/1e3+5, 5, "5 m/s",coordinates='data',color='k', alpha=0.8)
plt.xlabel('x (km)')
plt.ylabel('y (km)')
if db == 1:
    plt.title('{:s} at day = {:3.2f}'.format(titlename, times2D1[t1]))
    #plt.title('{:s} [{:s}] at day = {:3.2f}, W$_c$ = {:3.2f} m/s'.format(varname, units, times2D[t2], W_crit))
    #plt.title('mid-level {:s}, day = {:3.2f}'.format(varname, times2D[t2]))
else:
    plt.title('{:s} at day = {:3.2f}, block-averaging over ({:2.0f} km)$^2$'.format(titlename, times2D1[t1], db*(np.diff(x1)[0])/1e3))
    #plt.title('mid-level {:s}, day = {:3.2f}, block-averaging over ({:2.0f} km)$^2$'.format(varname, times2D[t2], db*(np.diff(x)[0])/1e3))
    #plt.title('{:s} [{:s}] at p = {:3.0f} hPa, at t = {:3.2f} days, block-averaging over ({:2.0f} km)$^2$'.format(varname, units, p[plevi], times2D[t2], db*(np.diff(x)[0])/1e3))
#cb = plt.colorbar(ticks=np.arange(10,90), 10)

#plt.legend(loc='best')
plt.savefig(fout + 'TILED{:s}map_day{:2.2f}db{:2.0f}.jpg'.format(varname, times2D1[t1], db))
#plt.savefig(fout + '{:s}map_p{:3.0f}_day{:2.1f}new.jpg'.format(varname, p[plevi], times2D[t2]))
plt.close()

plt.figure(2)
cs = plt.contourf(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, PW_field2new[:,:], PWvals, cmap=cm.RdYlBu_r, zorder=0)
ax=plt.gca()
#ticklabels = np.arange(0, 3500, 500)
ax.set_xticks(major_ticks2)
ax.set_yticks(major_ticks2)
minorlocator = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(minorlocator)
ax.yaxis.set_minor_locator(minorlocator)
ax.grid(which='major', color='k', zorder=10)
#ax.grid(which='minor', color='k', zorder=10)
#plt.axhline(1536, color='k', linewidth=0.5, alpha=0.5)
#plt.axvline(1536, color='k', linewidth=0.5, alpha=0.65)
#cs.collections[0].set_label('W$_c$ = {:3.2f} m/s'.format(W_crit))
#cs.collections[0].set_label('PW')
#plt.legend(loc='best')
#q = plt.quiver(xx[::8, ::8]/1e3, yy[::8, ::8]/1e3, U_tave[::8,::8], V_tave[::8,::8], scale=500, alpha=0.8, zorder=1)
if varname == 'BLQVADV':
    xxplot = xx[:-1, :-1]
    yyplot = yy[:-1, :-1]
else:
    xxplot = xx
    yyplot = yy
#plt.contourf(xxplot[::db, ::db]/1e3, yyplot[::db, ::db]/1e3, field_tave[:,:], vals, cmap=cm.RdYlBu_r, zorder=0)
cb = plt.colorbar()
cb.set_label('({0})'.format(units))
#p = plt.quiverkey(q, np.min(x)/1e3+30, np.max(y)/1e3+5, 5, "5 m/s",coordinates='data',color='k', alpha=0.8)
plt.xlabel('x (km)')
plt.ylabel('y (km)')
if db == 1:
    plt.title('{:s} at day = {:3.2f}'.format(titlename, times2D2[t2]))
    #plt.title('{:s} [{:s}] at day = {:3.2f}, W$_c$ = {:3.2f} m/s'.format(varname, units, times2D[t2], W_crit))
    #plt.title('mid-level {:s}, day = {:3.2f}'.format(varname, times2D[t2]))
else:
    plt.title('{:s} at day = {:3.2f}, block-averaging over ({:2.0f} km)$^2$'.format(titlename, times2D2[t2], db*(np.diff(x1)[0])/1e3))
    #plt.title('mid-level {:s}, day = {:3.2f}, block-averaging over ({:2.0f} km)$^2$'.format(varname, times2D[t2], db*(np.diff(x)[0])/1e3))
    #plt.title('{:s} [{:s}] at p = {:3.0f} hPa, at t = {:3.2f} days, block-averaging over ({:2.0f} km)$^2$'.format(varname, units, p[plevi], times2D[t2], db*(np.diff(x)[0])/1e3))
#cb = plt.colorbar(ticks=np.arange(10,90), 10)

#plt.legend(loc='best')
plt.savefig(fout + 'TILED{:s}map_day{:2.2f}db{:2.0f}.jpg'.format(varname, times2D2[t2], db))
#plt.savefig(fout + '{:s}map_p{:3.0f}_day{:2.1f}new.jpg'.format(varname, p[plevi], times2D[t2]))
plt.close()

plt.figure(3)
cs = plt.contourf(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, PW_field3[:,:], PWvals, cmap=cm.RdYlBu_r, zorder=0)
ax=plt.gca()
#ticklabels = np.arange(0, 3500, 500)
ax.set_xticks(major_ticks2)
ax.set_yticks(major_ticks2)
minorlocator = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(minorlocator)
ax.yaxis.set_minor_locator(minorlocator)
ax.grid(which='major', color='k', zorder=10)
#ax.grid(which='minor', color='k', zorder=10)
#plt.axhline(1536, color='k', linewidth=0.5, alpha=0.5)
#plt.axvline(1536, color='k', linewidth=0.5, alpha=0.65)
#cs.collections[0].set_label('W$_c$ = {:3.2f} m/s'.format(W_crit))
#cs.collections[0].set_label('PW')
#plt.legend(loc='best')
#q = plt.quiver(xx[::8, ::8]/1e3, yy[::8, ::8]/1e3, U_tave[::8,::8], V_tave[::8,::8], scale=500, alpha=0.8, zorder=1)
if varname == 'BLQVADV':
    xxplot = xx[:-1, :-1]
    yyplot = yy[:-1, :-1]
else:
    xxplot = xx
    yyplot = yy
#plt.contourf(xxplot[::db, ::db]/1e3, yyplot[::db, ::db]/1e3, field_tave[:,:], vals, cmap=cm.RdYlBu_r, zorder=0)
cb = plt.colorbar()
cb.set_label('({0})'.format(units))
#p = plt.quiverkey(q, np.min(x)/1e3+30, np.max(y)/1e3+5, 5, "5 m/s",coordinates='data',color='k', alpha=0.8)
plt.xlabel('x (km)')
plt.ylabel('y (km)')
if db == 1:
    plt.title('{:s} at day = {:3.2f}'.format(titlename, times2D3[t2]))
    #plt.title('{:s} [{:s}] at day = {:3.2f}, W$_c$ = {:3.2f} m/s'.format(varname, units, times2D[t2], W_crit))
    #plt.title('mid-level {:s}, day = {:3.2f}'.format(varname, times2D[t2]))
else:
    plt.title('{:s} at day = {:3.2f}, block-averaging over ({:2.0f} km)$^2$'.format(titlename, times2D3[t2], db*(np.diff(x1)[0])/1e3))
    #plt.title('mid-level {:s}, day = {:3.2f}, block-averaging over ({:2.0f} km)$^2$'.format(varname, times2D[t2], db*(np.diff(x)[0])/1e3))
    #plt.title('{:s} [{:s}] at p = {:3.0f} hPa, at t = {:3.2f} days, block-averaging over ({:2.0f} km)$^2$'.format(varname, units, p[plevi], times2D[t2], db*(np.diff(x)[0])/1e3))
#cb = plt.colorbar(ticks=np.arange(10,90), 10)

#plt.legend(loc='best')
plt.savefig(fout + 'TILED{:s}map_day{:2.2f}db{:2.0f}.jpg'.format(varname, times2D3[t2], db))
#plt.savefig(fout + '{:s}map_p{:3.0f}_day{:2.1f}new.jpg'.format(varname, p[plevi], times2D[t2]))
plt.close()
    
#plt.figure(2)
#plt.plot(times, alpha_t)
#plt.title('r$\overline{\alpha}$')
#plt.xlabel('time (days)')
#plt.ylabel('r$\overline{\alpha}$')
#plt.savefig(fout + '{:s}bar_day{:2.2f}db{:2.0f}.jpg'.format(varname, times2D[times[0]], times2D[times[-1]], db))








