from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from thermolib.wsat import wsat
from thermolib.constants import constants
from SAM_init_plot.block_fns import blockave3D
import SAM_init_plot.misc_fns 

c = constants()

matplotlib.rcParams.update({'font.size': 26})
matplotlib.rcParams.update({'figure.figsize': (20, 13)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 22})

plt.style.use('seaborn-white')

fpath2D =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpathSTAT = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'

foutdata = '/Users/cpatrizio/data/SST302/'

#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to140_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to130_64vert_ubarzero_MOISTDRYPROFS/'

fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_MOISTPROFS/'

nc_in2D1 = glob.glob(fpath2D + '*256x256*3000m*day230to250*302K.nc')[0]
nc_in3D1 = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]

nc_in2D2 = glob.glob(fpath2D + '*512x512*3000m*day180to195*302K.nc')[0]
nc_in3D2 = glob.glob(fpath3D + '*512x512*3000m*day180to195*302K.nc')[0]

nc_in2D3 = glob.glob(fpath2D + '*1024x1024*3000m*day170to180*302K.nc')[0]
nc_in3D3 = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]

domsizes = [768, 1536, 3072]
#nc_STATs = [nc_inSTAT1, nc_inSTAT2, nc_inSTAT3]
nc_2Ds = [nc_in2D1, nc_in2D2, nc_in2D3]
nc_3Ds = [nc_in3D1, nc_in3D2, nc_in3D3]

colors = ['k', 'r', 'g']

nave=10
ntave2D=24
ntave3D=4
#t2 = -35*ntave2D
#t3 = -35*ntave3D

t2s = [-1, -1, -1]
t3s = [-1, -1, -1]


for i, domsize in enumerate(domsizes):
    
    print 'domain size', domsize
    if domsize == 768:
        nave=10
    elif domsize == 1536:
        nave=10
    else:
        nave=5
        
    aveperiod2D = nave*ntave2D
    aveperiod3D = nave*ntave3D
    
    #nc_inSTAT = nc_STATs[i]
    nc_in2D = nc_2Ds[i]
    nc_in3D = nc_3Ds[i]
    t2=t2s[i]
    t3=t3s[i]
    
    print 'loading 2D variables'
    nc_data2D = Dataset(nc_in2D)
    varis2D = nc_data2D.variables
    
    PW = varis2D['PW'][t2-aveperiod2D:t2,:,:]
    
    print 'loading 3D variables'
    nc_data3D = Dataset(nc_in3D)
    varis3D = nc_data3D.variables
    
    t3D = varis3D['time'][t3-aveperiod3D:t3]
    t2D = varis2D['time'][t2-aveperiod2D:t2]
    x = varis3D['x'][:]
    y = varis3D['y'][:]
    z = varis3D['z'][:]
    p = varis3D['p'][:]
    p = p*1e2
    
    #calculate relative humidity 
    
    #averaging time period for 2D & 3D fields
    ntave2D=24
    ntave3D=4
    
    nt3D = t3D.size
    nt2D = t2D.size
    nz = z.size
    nx = x.size
    ny = y.size
    
    xx, yy = np.meshgrid(x, y)
    times = np.arange(t3D[0], np.max(t3D))
    
    QN = varis3D['QN'][t3-aveperiod3D:t3,:,:,:]
    QN_tave = np.mean(QN, axis=0)
    cldamnt = np.zeros(z.shape)
    
    p_t = 100
    z_t = z[p <= p_t*1e2][0]

    totpoints=nx*ny
    
    for zi, zlev in enumerate(z):
        QNz = QN_tave[zi,:,:]
        cldpoints = len(QNz[QNz > 0])
        cldamnt[zi] = cldpoints/(1.*totpoints)
        
    plt.figure(1)
    plt.plot(cldamnt, z/1e3, color=colors[i], label='{:d} km, day {:3.0f} to {:3.0f} average'.format(domsize, t3D[0], t3D[-1]))
    plt.xlabel('cloud fraction')
    plt.ylabel('z (km)')
    plt.ylim((0, z_t/1e3))
    plt.title(r'cloud fraction ($q_n$ > 0)')
    plt.savefig(fout + 'cldfrac_day250{:d}daynew.pdf'.format(nave))  
            

plt.legend(loc='best')
plt.savefig(fout + 'cldfrac_day250_{:d}daynew.pdf'.format(nave))    
plt.close()


