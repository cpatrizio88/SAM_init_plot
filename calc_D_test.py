from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from SAM_init_plot.block_fns import blockave2D, blockave3D
from thermolib.constants import constants
import matplotlib.colors as colors

c = constants()

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fpathSTAT = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to170_64vert_ubarzero_FREQDIST/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday90to160_64vert_ubarzero_FREQDIST/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STAT/'

#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]

#nc_in2D = glob.glob(fpath2D + '*1024x1024*3000m*day090to130*302K.nc')[0]


#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to140_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to130_64vert_ubarzero_MOISTDRYPROFS/'

fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_MOISTDRYPROFS/'

nc_inSTAT1 = glob.glob(fpathSTAT + '*256x256*3000m*250days*302K.nc')[0]
nc_in2D1 = glob.glob(fpath2D + '*256x256*3000m*day230to250*302K.nc')[0]
nc_in3D1 = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]


nc_inSTAT2 = glob.glob(fpathSTAT + '*512x512*3000m*180days*302K.nc')[0]
nc_in2D2 = glob.glob(fpath2D + '*512x512*3000m*day180to195*302K.nc')[0]
nc_in3D2 = glob.glob(fpath3D + '*512x512*3000m*day180to195*302K.nc')[0]

nc_inSTAT3 = glob.glob(fpathSTAT + '*1024x1024*3000m*180days*302K.nc')[0]
nc_in2D3 = glob.glob(fpath2D + '*1024x1024*3000m*day170to180*302K.nc')[0]
nc_in3D3 = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]


domsizes = [768, 1536, 3072]
nc_STATs = [nc_inSTAT1, nc_inSTAT2, nc_inSTAT3]
nc_2Ds = [nc_in2D1, nc_in2D2, nc_in2D3]
nc_3Ds = [nc_in3D1, nc_in3D2, nc_in3D3]

#nc_3Ds = [nc_in3D2]

pcolors = ['k', 'r', 'g']


tstarts = [0,90,90]

for i, nc_f3D in enumerate(nc_3Ds):
        
    print 'domsize', domsizes[i]
    
    print 'loading netCDF files'
    
    convfracs = np.array([])
    
    db=16
    W_c = 0.01
    
    times = np.array([])
    
    varname='W'

    #for nc_in3D in nc_f3D:
    #nc_data3D = Dataset(nc_in3D) 
        
    print 'loading', nc_f3D
    
    nc_data3D = Dataset(nc_f3D)
    varis3D = nc_data3D.variables
    
    nc_data2D = Dataset(nc_2Ds[i])#
    varis2D = nc_data2D.variables
    
    nc_dataSTAT = Dataset(nc_STATs[i])
    varisSTAT = nc_dataSTAT.variables
    
    P = varisSTAT['PREC'][:]
    
    rho_w = 1000
    #convert from mm/day to kg m^-2 s^-1
    P = P*(86400/1000)*rho_w
    
    x = varis3D['x'][:]
    y = varis3D['y'][:]
    z = varisSTAT['z'][:]
    p = varisSTAT['p'][:]
    p=p*1e2
    
    
    #nt3D = t3D.size
    nz = z.size
    nx = x.size
    ny = y.size
    
    nxprime = nx / db
    nyprime = ny / db
    
    #averaging time period for 2D & 3D fields
    ntave2D=24
    ntave3D=4

    t3D = varis3D['time'][:]
    nt2D = t3D.size
    field = varis3D[varname][:]
    
    P2D = varis2D['Prec'][:]
    
    #ntrunc = nt2D%ntave2D
    #ntrunc=0
    #field = field[ntrunc:,:,:]
    units = varis3D[varname].units
    #field_tmp = field.reshape(nt2D/ntave2D, ntave2D, nx, ny)
    #field_tave = np.mean(field_tmp, axis=1)
    nt = field.shape[0]
    convfrac = np.zeros(nt)
    #subsfrac = np.zeros(nt)
    for k, t in enumerate(t3D):
        k2 = (ntave2D/ntave3D)*k
        #Pfield_t = P2D[k2,:,:]
        field_t = field[k,:,:,:]
        #w = varis3D['W'][k,:,:,:]
        #T = varis3D['TABS'][k,:,:,:]
        #nt = T.shape[0]
        #p3D = np.zeros((ny, nx, nz))
        #p3D[:,:,:] = p
        #p3D = p3D[:,:,:,np.newaxis]
        #rho = p3D/(T*c.Rd)
        #massflux = np.multiply(w, rho)
        #massflux_blockave = blockave3D(massflux, db)
        if db == 1:
            field_blockave = field_t
            #P_blockave = Pfield_t
        else:
            field_blockave = blockave3D(field_t, db)
            #P_blockave = blockave2D(Pfield_t, db)
        field_vertave = np.mean(field_blockave, axis=0) 
        #xx, yy = np.meshgrid(x, y)
        #plt.figure(1)
        #plt.contour(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, field_vertave, [W_c], colos='k')
        #plt.contourf(xx[::db,::db]/1e3, yy[::db, ::db]/1e3, field_vertave, 50, cmap='RdYlBu_r')
        #plt.colorbar()
        #plt.xlabel('x (km)')
        #plt.ylabel('y (km')
        #plt.savefig(fout + '{:d}km_vertaveWmap_t{:3.0f}.pdf'.format(domsizes[i], t))
        #plt.close()
        #plt.figure(1)
        #plt.contour(xx[::db, ::db]/1e3, yy[::db, ::db]/1e3, field_vertave, [W_c], colos='k')
        #plt.contourf(xx[::db,::db]/1e3, yy[::db, ::db]/1e3, P_blockave, 50, cmap='GnBu_r')
        #plt.colorbar()
        #plt.xlabel('x (km)')
        #plt.ylabel('y (km')
        #plt.savefig(fout + '{:d}km_Pmap_t{:3.0f}.pdf'.format(domsizes[i], t))
        #plt.close()
        convfrac[k] = len(field_vertave[field_vertave[:,:] > W_c])/(1.*nxprime*nyprime)
        #subsfrac[i] = 1 - convfrac[i]
    
    #        buoyflux = np.vstack((buoyflux, buoyflux_temp)) if buoyflux.size else buoyflux_temp
    convfracs = np.concatenate((convfracs, convfrac)) if convfracs.size else convfrac
    times = np.concatenate((times, t3D)) if times.size else t3D
    
    #ntimes = convfracs.shape[0]
    
    print 'plotting'
    
    #calculate 2D (CRH rank, time) map of CRH frequency  ###
    #times = np.arange(tstarts[i], tstarts[i]+ntimes)

    w = 24*5
    
    convfracs_smooth = moving_average(convfracs, w)
    
    plt.figure(1)
    plt.plot(times, convfracs, color=pcolors[i], alpha=0.5, label='{:d} km'.format(domsizes[i]))
    plt.plot(times[w/2:-w/2+1], convfracs_smooth, color=pcolors[i], linewidth=1.5)
    #plt.plot(times, 1-convfracs, 'k--')
    plt.xlabel('time (days)')
    plt.ylabel('convective fractional area')
    plt.ylim((0, 0.10))
    if db == 1:
        plt.title('convective fractional area (vertically averaged W  > W$_c$ = {:3.2f} m/s)'.format(W_c))
    else:
        plt.title('convective fractional area (vertically averaged W > W$_c$ = {:3.2f} m/s), average over ({:2.0f} km)$^2$ blocks '.format(W_c, db*np.diff(x)[0]/1e3))
    plt.savefig(fout + 'fracconvNEW250days_db{:d}.pdf'.format(db))
  

plt.legend()
plt.savefig(fout + 'fracconv250days_db{:d}.pdf'.format(db))
plt.close()
        

    