# import matplotlib as mpl
# mpl.use('Agg')
from netCDF4 import Dataset
import site
site.addsitedir('/Users/cpatrizio/repos/thermolib/')
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from constants import constants

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

c=constants()

matplotlib.rcParams.update({'font.size': 28})
matplotlib.rcParams.update({'figure.figsize': (18, 10)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 22})

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
#foutSTAT = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_STAT/'
#foutSTAT = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggr130days_64vert_ubarzero_STAT/'
fout='/Volumes/GoogleDrive/My Drive/MS/MS figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STAT/'
fsave = '/Users/cpatrizio/data/SAM/'

# fpath2D = '/glade/scratch/patrizio/OUT_2D_nc/'
# fpath3D = '/glade/scratch/patrizio/OUT_3D_nc/'
# fout = '/glade/scratch/patrizio/OUT_STAT_FIGS/'


nc_in3D1 = glob.glob(fpath3D + '*256x256*day230to250*.nc')
nc_in3D2 = glob.glob(fpath3D + '*512x512*day180to195*.nc')
nc_in3D3 = glob.glob(fpath3D + '*1024x1024*day170to180*.nc')
nc_in3D4 = glob.glob(fpath3D + '*2048x2048*.nc')

nc_in2D1 = glob.glob(fpath2D + '*256x256*day230to250*.nc')
nc_in2D2 = glob.glob(fpath2D + '*512x512*day180to195*.nc')
nc_in2D3 = glob.glob(fpath2D + '*1024x1024*day170to180*.nc')
nc_in2D4 = glob.glob(fpath2D + '*2048x2048*.nc')

domsizes = [1536, 3072]
domsizes = [3072]
#domsizes = [768, 1536, 3072]
#domsizes = [768]
#domsizes = [768, 1536]
#domsizes = [3072]
#domsizes = [768, 1536, 3072]
#domsizes = [1536]
# domsizes = [768]
#domsizes = [1536]
#domsizes = [768, 1536, 3072, 6144]
#nc_3Ds = [nc_in3D2, nc_in3D3, nc_in3D4]
#nc_2Ds = [nc_in2D2, nc_in2D3, nc_in2D4]
nc_3Ds = [nc_in3D3]
nc_2Ds = [nc_in2D3]
# nc_3Ds = [nc_in3D2]
# nc_2Ds = [nc_in2D2]
#nc_3Ds = [nc_in3D1, nc_in3D2]
#nc_2Ds = [nc_in2D1, nc_in2D2]
#nc_3Ds = [nc_in3D1, nc_in3D2, nc_in3D3]
#nc_2Ds = [nc_in2D1, nc_in2D2, nc_in2D3]
#colors = ['r', 'g', 'm']
colors = ['g']
#colors = ['r', 'g']
#colors = ['k', 'r']

#domsize=768
#domsize=1536
#domsize=3072

print 'nc_in3D1', nc_in3D1
print 'nc_in3D2', nc_in3D2
print 'nc_in3D3', nc_in3D3
print 'nc_in3D4', nc_in3D4

print 'nc_in2D1', nc_in2D1
print 'nc_in2D2', nc_in2D2
print 'nc_in2D3', nc_in2D3
print 'nc_in2D4', nc_in2D4

#loop through domain cases
for j, domsize in enumerate(domsizes):
    
    print 'domsize', domsize
    nc_3Ddom = nc_3Ds[j]
    nc_2Ddom = nc_2Ds[j]
    
    #loop through each 3D field time series 
    for k, nc_3D in enumerate(nc_3Ddom):
        
        end = len(nc_3Ddom)
    
        nc_2D = nc_2Ddom[k]
        print '2D nc file', nc_2D
        print '3D nc file', nc_3D

        nc_data3D = Dataset(nc_3D)
        nc_data2D = Dataset(nc_2D)
        varis3D = nc_data3D.variables
        varis2D = nc_data2D.variables
        
        t3D = varis3D['time'][:]
        t2D = varis2D['time'][:]
        x = varis3D['x'][:]
        y = varis3D['y'][:]
        z = varis3D['z'][:]
        p = varis3D['p'][:]
        p = p*1e2
        
        dt = np.diff(t3D)[0]*(24*3600)
        delx = np.diff(x)[0]
        dely = np.diff(y)[0]
        
        nt3D = t3D.size
        nt2D = t2D.size
        nz = z.size
        nx = x.size
        ny = y.size
        
        qv = varis3D['QV'][:]
        qv=qv*1e-3
        qn = varis3D['QN'][:]
        qn=qn*1e-3
        T = varis3D['TABS'][:]
        
        P = varis2D['Prec'][:]
        LHF = varis2D['LHF'][:]
        SHF = varis2D['SHF'][:]
        
        qnet_TOA = varis2D['SWNT'][:] - varis2D['LWNT'][:]
        qnet_s = varis2D['SWNS'][:] - varis2D['LWNS'][:]
        
        Enet = (qnet_TOA-qnet_s) + LHF + SHF
        
        ntrunc = nt2D % 6

        #t2Di = t2D[0]
        #t3Di = t3D[0]

        #ix =np.argmin(np.abs(t2D == t3Di))
        
        Enet = Enet[ntrunc:,:,:]
        P = P[ntrunc:,:,:]
        qnet_TOA = qnet_TOA[ntrunc:,:,:]
        qnet_s = qnet_s[ntrunc:,:,:]
        #P = P[::6,:,:]
        P = P.reshape(nt2D/6, 6, nx, ny)
        P = np.mean(P, axis=1)

        Enet = Enet.reshape(nt2D/6, 6, nx, ny)
        Enet = np.mean(Enet, axis=1)
        qnet_TOA = qnet_TOA.reshape(nt2D/6, 6, nx, ny)
        qnet_TOA = np.mean(qnet_TOA, axis=1)
        qnet_s = qnet_s.reshape(nt2D/6, 6, nx, ny)
        qnet_s = np.mean(qnet_s, axis=1)
        #Enet = Enet[::6,:,:]
        #qnet_TOA = qnet_TOA[::6,:,:]
        #qnet_s = qnet_s[::6,:,:]
        print 'length t2D', len(t2D)
        print 'length t3D', len(t3D)  
        print 'length of 2D array after reshaping', qnet_s.shape[0]

        z3D = np.zeros((nz, nx, ny)).T
        z3D[:,:,:] = z
        z3D = z3D.T
        
       
        
        rhoz = p/(c.Rd*np.mean(np.mean(T, axis=3), axis=2))
            
        rhoz_tave = np.mean(rhoz, axis=0)
        rhoz3D = np.zeros((ny, nx, nz))
        rhoz3D[:,:,:] = rhoz_tave
        rhoz3D = rhoz3D.T
        
        delz = np.diff(z)
        delz3D = np.zeros((ny, nx, nz-1))
        delz3D[:,:,:] = delz
        delz3D = delz3D.T

        delp = np.diff(p)
        delp3D = np.zeros((ny, nx, nz-1))
        delp3D[:,:,:] = delp
        delp3D = delp3D.T
        
        print 'max qn (kg/kg)', np.ma.max(qn)
        
        #TEST LIQUID ICE STATIC ENERGY 
        qi = np.ma.masked_array(qn, mask = T > 273)
        np.ma.set_fill_value(qi, 0)
        qi = qi.data
        print 'max qi (kg/kg)', np.ma.max(qi)
        #qi = np.zeros(qn.shape) 
        #TEST MOIST STATIC ENERGY
        #qi = np.zeros(qn.shape)
        Plen = P.shape[0] 
        dhhatdt = np.zeros((len(t3D)-2, nx, ny))
        hhat = np.zeros((len(t3D)-2, nx, ny))

        for i, t in enumerate(np.arange(1,len(t3D)-1)):
            print 'day', t3D[t]
            h_t1 = c.cpd*T[t-1,:] + c.g*z3D + c.lv0*qv[t-1,:] - c.lf*qi[t-1,:]
            h_t2 = c.cpd*T[t+1,:] + c.g*z3D + c.lv0*qv[t+1,:] - c.lf*qi[t+1,:]
            
            h_t1hat = np.sum(h_t1[:-1,:]*rhoz3D[:-1,:]*delz3D, axis=0)
            h_t2hat = np.sum(h_t2[:-1,:]*rhoz3D[:-1,:]*delz3D, axis=0)
            #h_t1hat = np.nansum(h_t1[:-1,:]*delp3D, axis=0)
            #h_t2hat = np.nansum(h_t2[:-1,:]*delp3D, axis=0)
            hhat[i, :, :] = (h_t2hat + h_t1hat)/2. 
            dhhatdt[i, :, :] = (h_t2hat - h_t1hat)/(2*dt)
        
        P_c = 0.0
    
        conv_points =  P[1:-1,:]  > P_c
        dry_points = P[1:-1,:] <= P_c
        
        hhatbar = np.mean(np.mean(hhat, axis=2), axis=1)
        dhhatdtbar = np.mean(np.mean(dhhatdt, axis=2), axis=1)
        Enetbar = np.mean(np.mean(Enet[1:-1,:], axis=2), axis=1)
        #A_c = np.sum(conv_points)*delx*dely
        #A_d = np.sum(dry_points)*delx*dely
        A_c = np.sum(np.sum(conv_points, axis=1),axis=1)*delx*dely 
        
        dhhatdt_d = np.ma.masked_array(dhhatdt, mask = conv_points)
        Enet_d = np.ma.masked_array(Enet[1:-1,:], mask = conv_points)
        
        dhhatdt_c = np.ma.masked_array(dhhatdt, mask = dry_points)
        Enet_c = np.ma.masked_array(Enet[1:-1,:], mask = dry_points)
       
        qnet_s_c = np.ma.masked_array(qnet_s[1:-1,:], mask = dry_points)
        qnet_TOA_c = np.ma.masked_array(qnet_TOA[1:-1,:], mask = dry_points)

        qnet_s_d = np.ma.masked_array(qnet_s[1:-1,:], mask = conv_points)
        qnet_TOA_d = np.ma.masked_array(qnet_TOA[1:-1,:], mask = conv_points)

        qnet_s_d = np.nanmean(np.nanmean(qnet_s_d, axis=2), axis=1)
        qnet_TOA_d = np.nanmean(np.nanmean(qnet_TOA_d, axis=2), axis=1)

        qnet_s_c = np.nanmean(np.nanmean(qnet_s_c, axis=2), axis=1)
        qnet_TOA_c = np.nanmean(np.nanmean(qnet_TOA_c, axis=2), axis=1)

        dhhatdt_cbar = np.nanmean(np.nanmean(dhhatdt_c, axis=2), axis=1)
        dhhatdt_dbar = np.nanmean(np.nanmean(dhhatdt_d, axis=2), axis=1)

        divh_d = np.nansum(np.nansum(Enet_d - dhhatdt_d, axis=2), axis=1)*delx*dely
        divh_c = np.nansum(np.nansum(Enet_c - dhhatdt_c, axis=2), axis=1)*delx*dely
       
        #divh_d = np.nanmean(np.nanmean(Enet_d - dhhatdt_d, axis=2), axis=1)
        #divh_c = np.nanmean(np.nanmean(Enet_c - dhhatdt_c, axis=2), axis=1)
        hhatbar = np.ma.filled(hhatbar, np.nan)
        Enetbar = np.ma.filled(Enetbar, np.nan)
        dhhatdtbar = np.ma.filled(dhhatdtbar, np.nan) 
        divh_c = np.ma.filled(divh_c, np.nan)
        divh_d = np.ma.filled(divh_d, np.nan)
        qnet_s_c = np.ma.filled(qnet_s_c, np.nan)
        qnet_TOA_c = np.ma.filled(qnet_TOA_c, np.nan)
        qnet_s_d = np.ma.filled(qnet_s_d, np.nan)
        qnet_TOA_d = np.ma.filled(qnet_TOA_d, np.nan)
        dhhatdt_cbar = np.ma.filled(dhhatdt_cbar, np.nan)
        dhhatdt_dbar = np.ma.filled(dhhatdt_dbar, np.nan)

        N = 5*4
        
        smooth_divh_c = running_mean(divh_c, N)
        t3D = t3D[1:-1]
        
        fname = '{:d}km_div_h_c_day{:3.0f}to{:3.0f}_{:2.1f}mm.npz'.format(domsize, t3D[0], t3D[-1], P_c)
        np.savez(fsave + fname, dhhatdtbar=dhhatdtbar, hhatbar=hhatbar, Enetbar=Enetbar, divh=divh_c, qnets=qnet_s_c, qnetTOA=qnet_TOA_c, dhhatdt=dhhatdt_cbar, time=t3D, A=A_c)

        fname2 = '{:d}km_div_h_d_day{:3.0f}to{:3.0f}_{:2.1f}mm.npz'.format(domsize, t3D[0], t3D[-1], P_c)
        np.savez(fsave + fname2, divh=divh_d, qnets=qnet_s_d, qnetTOA=qnet_TOA_d, dhhatdt=dhhatdt_dbar, time=t3D)
        
        

        plt.figure(1)
        #plt.plot(t3D[N/2-1:-N/2], -smooth_divh_c, color=colors[j], label='{:d} km'.format(domsize))
        #plt.plot(t3D, -divh_c, color=colors[j], linewidth=2, alpha=0.5)
        if k == end:
            plt.plot(t3D[N/2-1:-N/2], smooth_divh_c, color=colors[j], label='{:d} km'.format(domsize))
        else:
            plt.plot(t3D[N/2-1:-N/2], smooth_divh_c, color=colors[j])
        plt.plot(t3D, divh_c, color=colors[j], linewidth=2, alpha=0.5)
        plt.axhline(0, color='b')
        #plt.ylim(-50,200)
        plt.xlabel('time (days)')
        plt.title(r'$M$')
        plt.ylabel('total MSE export from convective region (J/s)')
        plt.savefig(fout + 'MSEexport_Pc{:1.1f}.pdf'.format(P_c))
        
        divh_c_norm = divh_c/A_c
        
        smooth_divh_c_norm = running_mean(divh_c_norm, N)
        
        plt.figure(2)
        if k == end:
            plt.plot(t3D[N/2-1:-N/2], smooth_divh_c_norm, color=colors[j], label='{:d} km'.format(domsize))
        else:
            plt.plot(t3D[N/2-1:-N/2], smooth_divh_c_norm, color=colors[j])
        plt.plot(t3D, divh_c_norm, color=colors[j], linewidth=2, alpha=0.5)
        plt.axhline(0, color='b')
        plt.xlabel('time (days)')
        plt.title(r'$M$')
        plt.ylabel(r'$M$ (W m$^{-2}$)')
        plt.ylim(-50,200)
        plt.savefig(fout + 'MSEexportnorm_Pc{:1.1f}.pdf'.format(P_c))
        
        print 'domsize', domsize
        print '20-day mean MSE export (W/m^2):', np.mean(divh_c_norm[-4*20:])

plt.figure(1)     
plt.legend()
plt.savefig(fout + 'MSEexport_Pc{:1.1f}.pdf'.format(P_c))
plt.close()

plt.figure(2)
plt.legend()
plt.savefig(fout + 'MSEexportnorm_Pc{:1.1f}.pdf'.format(P_c))
plt.close()

plt.close('all')
    
