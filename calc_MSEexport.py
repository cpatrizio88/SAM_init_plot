# -*- coding: utf-8 -*-
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from thermolib.constants import constants

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

nc_in3D1 = glob.glob(fpath3D + '*256x256*302K.nc')
nc_in3D2 = glob.glob(fpath3D + '*512x512*3000m*day180to195*302K.nc')
nc_in3D3 = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]

nc_in2D1 = glob.glob(fpath2D + '*256x256*302K.nc')
nc_in2D2 = glob.glob(fpath2D + '*day180to195*302K.nc')
nc_in2D3 = glob.glob(fpath2D + '*day170to180*302K.nc')[0]

domsizes = [768, 1536]
nc_3Ds = [nc_in3D1, nc_in3D2, nc_in3D3]
nc_2Ds = [nc_in2D1, nc_in2D2, nc_in2D3]
colors = ['k', 'r', 'g']

#domsize=768
#domsize=1536
#domsize=3072

#loop through domain cases
for j, domsize in enumerate(domsizes):
    
    print 'domsize', domsize
    nc_3Ddom = nc_3Ds[j]
    nc_2Ddom = nc_2Ds[j]
    
    #loop through each 3D field time series 
    for k, nc_3D in enumerate(nc_3Ddom):
        
        end = len(nc_3Ddom)
    
        nc_2D = nc_2Ddom[k]
        print 'nc file', nc_2D

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
        qnet_s = varis2D['LWNS'][:] - varis2D['SWNS'][:]
        
        Enet_s = qnet_s + LHF + SHF
        
        ntrunc = nt2D % 6
        
        Enet_s = Enet_s[ntrunc:,:,:]
        P = P[ntrunc:,:,:]
        P = P.reshape(6, nt3D, nx, ny)
        P = np.mean(P, axis=0)
        
        Enet_s = Enet_s.reshape(6, nt3D, nx, ny)
        Enet_s = np.mean(Enet_s, axis=0)
        
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
        
        #TEST LIQUID ICE STATIC ENERGY 
        qi = np.ma.masked_array(qn, mask = T > 273)
        
        #TEST MOIST STATIC ENERGY
        #qi = np.zeros(qn.shape)
        
        dhhatdt = np.zeros((nt3D-1, nx, ny))
        
        for i, t in enumerate(t3D[:-1]):
            print 'day', t
            h_t1 = c.cpd*T[i,:] + c.g*z3D + c.lv0*qv[i,:] - c.lf*qi[i,:]
            h_t2 = c.cpd*T[i+1,:] + c.g*z3D + c.lv0*qv[i+1,:] - c.lf*qi[i+1,:]
            
            h_t1hat = np.nansum(h_t1[:-1,:]*rhoz3D[:-1,:]*delz3D, axis=0)
            h_t2hat = np.nansum(h_t2[:-1,:]*rhoz3D[:-1,:]*delz3D, axis=0)
            
            dhhatdt[i, :, :] = (h_t2hat - h_t1hat)/dt
            
        conv_points =  P[:-1,:] > 1
        dry_points = P[:-1,:] < 1
        
        dhhatdt_d = np.ma.masked_array(dhhatdt, mask = conv_points)
        Enet_s_d = np.ma.masked_array(Enet_s[:-1,:], mask = conv_points)
        
        dhhatdt_c = np.ma.masked_array(dhhatdt, mask = dry_points)
        Enet_s_c = np.ma.masked_array(Enet_s[:-1,:], mask = dry_points)
        
        divh_d = np.nansum(np.nansum(Enet_s_d*delx*dely - dhhatdt_d, axis=2), axis=1)
        divh_c = np.nansum(np.nansum(Enet_s_c*delx*dely - dhhatdt_c, axis=2), axis=1)
        
        N = 5*4
        
        smooth_divh_c = running_mean(divh_c, N)
        t3D = t3D[:-1]
        
        
        #plt.plot(t3D[N/2-1:-N/2], -smooth_divh, color=colors[j], label='{:d} km'.format(domsize))
        #plt.plot(t3D, -divh, color=colors[j], linewidth=2, alpha=0.5)
        if k == end:
            plt.plot(t3D[N/2-1:-N/2], smooth_divh_c, color=colors[j], label='{:d} km'.format(domsize))
        else:
            plt.plot(t3D[N/2-1:-N/2], smooth_divh_c, color=colors[j])
        plt.plot(t3D, divh_c, color=colors[j], linewidth=2, alpha=0.5)
        plt.axhline(0, color='b')
        plt.xlabel('time (days)')
        plt.ylabel('MSE export from convective region (J/s)')
        plt.savefig(fout + 'MSEexport.pdf')
    
plt.legend()
plt.savefig(fout + 'MSEexport.pdf')
plt.close()
    
