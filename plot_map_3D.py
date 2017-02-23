from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker
import matplotlib
from SAM_init_plot.block_fns import blockave3D, blockave2D
from thermolib.constants import constants
from thermolib.thermo import thermo
from thermolib.findTmoist import findTmoist
from thermolib.wsat import wsat

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
#nc_in = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
 
nc_STAT = glob.glob(fpathSTAT + '*512x512*3000m*180days*302K.nc')[0]
nc_in = glob.glob(fpath2D + '*512x512*3000m*day090to140*302K.nc')[0]
nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day090to130*302K.nc')[0]

#nc_STAT = glob.glob(fpathSTAT + '*1024x1024*3000m*170days*302K.nc')[0]
#nc_in = glob.glob(fpath2D + '*1024x1024*3000m*day090to130*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]

#domsize=768
domsize=1536
#domsize=3072

nc_data= Dataset(nc_in)
nc_dataSTAT = Dataset(nc_STAT)
nc_data3D = Dataset(nc_in3D)
varis = nc_data.variables
varis3D = nc_data3D.variables
varisSTAT = nc_dataSTAT.variables
times2D = varis['time'][:]

#times3D = varis3D['time'][:]

x = varis['x'][:]
y = varis['y'][:]

fac=8

nave2D = 24/fac
nave3D = 4

times = np.arange(0*fac, 40*fac)

db=1


xx, yy = np.meshgrid(x, y)
#LWNT: OLR
#Prec: surface precip
varname = 'W500'
vari = varis[varname]
#field = vari[:]
#field_tave = blockave3D(field, db)


p = varis3D['p'][:]
z = varis3D['z'][:]
varname = 'W'#
plev = 500
plevi = np.where(p > plev)[0][-1]
wfield_p = varis3D['W'][:,plevi,:,:]
#dbw=16
block_field = blockave3D(wfield_p, db)
field = block_field


nx=x.size
ny=y.size
nz=z.size


#tstep=int(4)
#tstep=int(12)
#tstart=4*90
#tstart = 24*40
#tstart=nave3D


W_crit = 0.5
Wedge_u = 3
Wedge_l = 0.5
vals=[Wedge_l*W_crit, 0.01, Wedge_u*W_crit]
vals=20

#CAPE COMPUTATION
varname = 'CAPE'
#TABS = varis3D['TABS'][:,:,:,:]
#QV = varis3D['QV'][:,:,:,:]
#QV = QV*1e-3
#
#TV = TABS*(1+0.61*QV)

#varname ='URADIAL'
#varname = 'Nsquared'
if varname == 'CAPE':
      units = r'J/kg'
elif varname == 'N':
      units = r's$^{{-2}}$'
elif varname == 'URADIAL':
      U_rfname = glob.glob(foutdata + '{:d}km_URADIAL_*{:3.0f}*'.format(domsize, times3D[-1]))[0] 
      U_r = np.load(U_rfname) 
      units = 'm/s'
else:
      units = vari.units.strip()

#delz = np.diff(z)

#vals = np.linspace(-2, 5, 50)

#PW values
#vals = np.arange(0, 88, 4)

#CWP values
#vals = np.arange(0, 1.3, 0.02)

#N^2 values
#vals = np.linspace(2.3e-4, 2.6e-4, 50)

#vals = np.linspace(-1, 1, 50)

#ZC vals
#vals = np.arange(0, 15, 0.2)


#CAPE vals?
vals = np.arange(0, 500, 10)




#for ii, i in enumerate(np.arange(tstart, len(times3D), tstep)):
#EDIT THIS TO PLOT PROJECTION OF 3D FIELD
for i in times:
    
    t2=i*nave2D
    t3=i*nave3D
    
    if varname == 'CAPE':
        
        #need to calculate z_n and z_f
     
        p_s = p[0]
        T_s = varisSTAT['TABS'][t2-nave2D:t2,0]
        q_sat = wsat(T_s, p_s*1e2)
        
        Tenv = varis3D['TABS'][t3-nave3D:t3,:,:,:]
        
        QV = varis3D['QV'][t3-nave3D:t3,:,:,:]
        QV = QV*1e-3
        
        Tenv_ave = np.mean(np.mean(Tenv, axis=2), axis=1)
        
        
        thetae0 = thermo.theta_e(T_s, p_s*1e2, q_sat, 0) #theta_e in moist region 
        
        Tadiabat = findTmoist(thetae0, p*1e2)
        
        qvadiabat = []
        
        for i, plev in enumerate(p*1e2):
            Tadb = Tadiabat[i]
            qvadiabat = qvadiabat + [wsat(Tadb, plev)]
            
        qvadiabat = np.array(qvadiabat)
            
        TVenv = Tenv*(1+0.61*QV)
        
        TVenv_ave = np.mean(np.mean(TVenv, axis=2), axis=1)
            
        TVadiabat = Tadiabat*(1+0.61*qvadiabat)
        
        
        TVadiabat3D = np.zeros((ny, nx, nz))
        TVadiabat3D[:,:,:] = TVadiabat
        TVadiabat3D = TVadiabat3D.T
        
        Tpert = TVadiabat3D - TVenv
        
        Tpert_ave = np.mean(np.mean(Tpert, axis=2), axis=1)
        
        znb_i = np.where(Tpert_ave[::-1] > 0)[0][0]
        z_nb = z[::-1][znb_i]
        znb_i = np.where(z > z_nb)[0][0]
        zL_i = np.where(Tpert_ave > 0)[0][0]
        
        zL_i = 0
        
        delz = np.diff(z)
        delz3D = np.zeros((ny, nx, nz-1))
        delz3D[:,:,:] = delz
        delz3D = delz3D.T
                
        CAPE = c.g*np.sum(np.multiply(Tpert[zL_i:znb_i,:,:]/TVenv[zL_i:znb_i,:,:], delz3D[zL_i:znb_i,:,:]), axis=0)
    if varname == 'Nsquared':
        p0=1000
        T = varis3D['TABS'][t3-nave3D:t3,:,:,:]
        T_tave = np.mean(T, axis=0)
        p3D = np.zeros((nz, nx, ny))
        p3D= p3D.T
        p3D[:,:,:] = p
        p3D = p3D.T
        #p3D = p3D.T
        #p4D = np.tile(p3D[:,:,:, np.newaxis], T.shape[0])
        #p4D = p4D.T
        THETA_tave = T_tave*(p0/p3D)**(c.Rd/c.cp)
        diffTHETA = THETA_tave[1:,:,:] - THETA_tave[:-1,:,:]
        diffz = np.diff(z)
        diffz3D = np.zeros((nz-1, nx, ny))
        diffz3D = diffz3D.T
        diffz3D[:,:,:] = diffz
        diffz3D = diffz3D.T
        THETA_mid = (THETA_tave[1:,:,:] + THETA_tave[:-1,:,:])/2.
        N_tave = (c.g*diffTHETA)/np.multiply(diffz3D, THETA_mid)
        N_tave_vertave = np.mean(N_tave, axis=0)
        field_tave = blockave2D(N_tave_vertave, db)
        
    if varname == 'URADIAL':
        U_r_p = U_r[plevi,:,:]
        field_tave = blockave2D(U_r_p, db)


    #w_tave = np.mean(wfield_p[t3-nave3D:t3,:,:], axis=0)
    #w_tave = blockave2D(w_tave, db)
    
    #wvals = [W_crit]
    
    

    
    #field_tave = np.mean(field[t2-nave2D:t2,:,:], axis=0)
    field_tave = blockave2D(field_tave, db)
    #varname_u = 'U200'
    #U = varis[varname_u][t2-nave2D:t2,:,:]
    #U_tave = np.mean(U, axis=0)
    #V = varis['V200'][i-nave2D:i,:,:]
    #V_tave = np.mean(V, axis=0)

    print 'day', times2D[t2]
    #print 'day', times3D[i]
    plt.figure()
    #minvar= np.min(field_tave[~np.isnan(field_tave)])
    #maxvar= np.max(field_tave[~np.isnan(field_tave)])
    #plt.contour(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, field_tave[:,:], vals, colors='k', linewidth=0.5, alpha=0.5)
    #plt.contour(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, w_tave[:,:], wvals, colors=('k',), alpha=0.8, linewidth=1, zorder=1, label = 'w at {:2.0f} hPa'.format(p[plevi]))
    #plt.legend(loc='best')
    #q = plt.quiver(xx[::8, ::8]/1e3, yy[::8, ::8]/1e3, U_tave[::8,::8], V_tave[::8,::8], scale=500, alpha=0.8, zorder=1)
    plt.contourf(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, field_tave[:,:], vals, cmap=cm.RdYlBu_r, zorder=0)
    cb = plt.colorbar()
    cb.set_label('({0})'.format(vari.units.strip()))
    #p = plt.quiverkey(q, np.min(x)/1e3+30, np.max(y)/1e3+5, 5, "5 m/s",coordinates='data',color='k', alpha=0.8)
    plt.xlabel('x (km)')
    plt.ylabel('y (km)')
    if db == 1:
        #plt.title('{:s} [{:s}] at day = {:3.2f}, W$_c$ = {:3.2f} m/s'.format(varname, units, times2D[t2], W_crit))
        plt.title('{:s} [{:s}] at day = {:3.2f}'.format(varname, units, times2D[t2]))
    else:
        plt.title('{:s} (m/s) and {:s} [{:s}] at t = {:3.2f} days, block-averaging over ({:2.0f} km)$^2$'.format(varname, units, times2D[i], db*(np.diff(x)[0])/1e3))
        #plt.title('{:s} [{:s}] at p = {:3.0f} hPa, at t = {:3.2f} days, block-averaging over ({:2.0f} km)$^2$'.format(varname, units, p[plevi], times2D[t2], db*(np.diff(x)[0])/1e3))
    #cb = plt.colorbar(ticks=np.arange(10,90), 10)
    #plt.legend(loc='best')
    #plt.savefig(fout + '{:s}map_day{:2.2f}db{:2.0f}.jpg'.format(varname, times2D[t2], db))
    plt.savefig(fout + '{:s}map_p{:3.0f}_day{:2.1f}new.jpg'.format(varname, p[plevi], times2D[t2]))
    plt.close()








