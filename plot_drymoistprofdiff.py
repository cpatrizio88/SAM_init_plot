from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from thermolib.wsat import wsat
from thermolib.constants import constants
import SAM_init_plot.misc_fns 

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (20, 13)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 16})

plt.style.use('seaborn-white')

fpath2D =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpathSTAT = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'

foutdata = '/Users/cpatrizio/data/SST302/'

#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to140_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to130_64vert_ubarzero_MOISTDRYPROFS/'

fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_MOISTDRYPROFS/'

nc_inSTAT1 = glob.glob(fpathSTAT + '*256x256*3000m*130days*302K.nc')[0]
nc_in2D1 = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]
nc_in3D1 = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]


nc_inSTAT2 = glob.glob(fpathSTAT + '*512x512*3000m*130days*302K.nc')[0]
nc_in2D2 = glob.glob(fpath2D + '*512x512*3000m*day90to130*302K.nc')[0]
nc_in3D2 = glob.glob(fpath3D + '*512x512*3000m*day90to130*302K.nc')[0]

nc_inSTAT3 = glob.glob(fpathSTAT + '*1024x1024*3000m*150days*302K.nc')[0]
nc_in2D3 = glob.glob(fpath2D + '*1024x1024*3000m*day140to150*302K.nc')[0]
nc_in3D3 = glob.glob(fpath3D + '*1024x1024*3000m*day140to150*302K.nc')[0]


domsizes = [768, 1536, 3072]
nc_STATs = [nc_inSTAT1, nc_inSTAT2, nc_inSTAT3]
nc_2Ds = [nc_in2D1, nc_in2D2, nc_in2D3]
nc_3Ds = [nc_in3D1, nc_in3D2, nc_in3D3]

colors = ['k', 'r', 'g']

nave=5
ntave2D=24
ntave3D=4
t2 = -35*ntave2D
t3 = -35*ntave3D

t2s = [-35*ntave2D, -1, -1]
t3s = [-35*ntave3D, -1, -1]

aveperiod2D = nave*ntave2D
aveperiod3D = nave*ntave3D

#varnames = ['QRAD', 'W', 'U', 'QN', 'QV']
#varnames=['QV', 'TABS']
#varnames=['QRAD', 'W', 'QN']
#varnames=['QV']
varnames=['TABS']

#MAYBE CALCULATE THETA, AND THEN STABILITY?

labels=[]
hs=[]
    
for i, domsize in enumerate(domsizes):
    
    print 'domain size', domsize
    
    nc_inSTAT = nc_STATs[i]
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
    
    print 'loading STAT variables'
    nc_dataSTAT = Dataset(nc_inSTAT)
    varisSTAT = nc_dataSTAT.variables
    
    theta = varisSTAT['THETA'][t2-aveperiod2D:t2,:]
    RELH = varisSTAT['RELH'][t2-aveperiod2D:t2,:]
    RELH_tave = np.mean(RELH, axis=0)
    theta_tave = np.mean(theta, axis=0)
    
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
    
    
    #EDIT HIS TO CHANGE VARIABLE FOR MOIST REGION THRESHOLD (2D threshold)
    mvarname = 'PW'
    mvar = varis2D[mvarname]
    mfield = mvar[t2-aveperiod2D:t2,:,:]
    mfield_t= np.mean(mfield, axis=0)

    #EDIT THIS TO CHANGE HOW MANY STD ABOVE MEAN FOR MOIST REGION THRESHOLD
    a=1.5
    
    mfieldcrit = np.mean(mfield) + a*np.std(mfield)
    moist_points = mfield_t >= mfieldcrit
    dry_points = mfield_t < mfieldcrit 
    
    #UNCOMMENT TO USE W > W_crit AS CONVECTIVE REGION THRESHOLD
    #W=varis3D['W'][t3-aveperiod3D:t3,:,:]
    #W_tave = np.mean(W, axis=0)
    #W_crit = 0.05
    #moist_points = W_tave >= W_crit
    #dry_points = W_tave < W_crit
    
    #find BL to get w_BL and w_m
    dthetadz = (theta_tave[1:]-theta_tave[:-1])/np.diff(z)
    
    dthetadz_crit = 1e-3
    
    BLi = np.where(dthetadz > dthetadz_crit)[0][0]
    
    p_BL = p[BLi]
    z_BL = z[BLi]
    zBL = z[:BLi]
    p_t = 150
    z_t = z[p <= p_t*1e2][0]
    
    
    print 'calculating dry region and convective region average profiles'
    print 'domain length = {:3.0f}, day {:3.0f} to day {:3.0f}'.format(domsize, t3D[0], t3D[-1])
    
    #plot radial velocity in moist region/convective region?
    
    #U_rfname = glob.glob(foutdata + '*{:d}*'.format(domsize))[0]
    
    #U_r = np.load(foutdata + U_rfname)
    
    for k, varname in enumerate(varnames):
        print 'loading', varname
        vari = varis3D[varname]
        field = varis3D[varname][t3-aveperiod3D:t3,:,:,:]
        field_tave = np.mean(field, axis=0)
        #calculate relative humidity 
        #if varname == 'TABS':
        #    varname = 'RH'
        #    QV = varis3D['QV'][t3-aveperiod3D:t3,:,:,:]
        #    QV_tave = np.mean(QV, axis=0)
        #    RH_tave = np.zeros(QV_tave.shape)
        #    for pi, plev in enumerate(p):
        #        wsat_tave = 1000*wsat(field_tave[pi,:,:], plev) #convert to g/kg
        #        RH_tave[pi,:,:] = 100*(QV_tave[pi,:,:]/wsat_tave) #convert to percent
        #    field_tave = RH_tave
        
        moist_field = field_tave[:,moist_points]
        dry_field = field_tave[:,dry_points]
        
        #moist_prof=[]
        #dry_prof=[]
        
        #print 'calculating convective/dry profile'
        #for zi, zlev in enumerate(z):
        #    moist_points = W_tave[zi,:,:] >= W_crit
        #    dry_points = W_tave[zi,:,:] < W_crit
        #    moist_bar = np.mean(field_tave[zi, moist_points])
        #    dry_bar = np.mean(field_tave[zi, dry_points])
        #    moist_prof= moist_prof + [moist_bar]
        #    dry_prof = dry_prof + [dry_bar]
        #    
        #moist_prof = np.array(moist_prof)
        #dry_prof = np.array(dry_prof)
            
        
        moist_prof = np.mean(moist_field, axis=1)
        dry_prof = np.mean(dry_field, axis=1)
        
        if varname == 'RH':
            units='%'
        else:
            units = vari.units.strip()
            
        interior = np.bitwise_and(z > z_BL, z < z_t)
        zint = z[interior]
        moist_prof = moist_prof[interior] 
        dry_prof = dry_prof[interior]
        
        moist_profave = np.mean(moist_prof)
        dry_profave = np.mean(dry_prof)
        
        fig = plt.figure(k)
        ax1 = fig.add_subplot(2,1,1)
        #plt.suptitle(r'{:s} ({:s}), convective region {:s} threshold = {:2.1f} {:s}'.format(varname, units, mvarname, mfieldcrit, mvar.units.strip()))
        h1, = ax1.plot(moist_prof, zint/1e3, color=colors[i])
        ax1.set_title('{:s}, convection region, using {:s} threshold'.format(varname, mvarname))
        #axarr[0,].set_xlabel('{:s} ({:s})'.format(varname, vari.units.strip()))
        ax1.set_ylabel('z (km)')
        ax1.set_ylim(z_BL/1e3, z_t/1e3)
        #if (varname != 'QV'):
        #    axarr[0,].axvline(moist_profave, color='k', alpha=0.5, label='vertical average')
        ax2 = fig.add_subplot(2,1,2)
        h2, = ax2.plot(dry_prof, zint/1e3, color=colors[i])
        ax2.set_title('{:s}, dry region'.format(varname))
        ax2.set_xlabel('{:s} ({:s})'.format(varname, vari.units.strip()))
        ax2.set_ylabel('z (km)')
        ax2.set_ylim(z_BL/1e3, z_t/1e3)
        #axarr[1,].axvline(dry_profave, color='k', alpha=0.5, label='vertical average')
        label1=r'day {:3.0f} to {:3.0f} average, domain size = {:d} km$^2$'.format(t3D[0], t3D[-1], domsize)
        #fig.legend(hs, labels, 'best')       
        fig.savefig(fout + '{:s}dryconvprof_{:s}thresh.pdf'.format(varname, mvarname))
    labels = labels + [label1]
    hs = hs + [h1]


for k, varname in enumerate(varnames):
    fig=plt.figure(k)
    plt.figlegend(hs, labels, 'best')
    fig.savefig(fout + '{:s}dryconvprof_{:s}thresh.pdf'.format(varname, mvarname))
    plt.close(k)
    
        
    ##simple model results   
    #w_BL = -.0066
    #w_m = 0.0334
    #
    #W_m = 3
    #W_d = -4.5
    #
    #q_FA = 0.000175*1e3 #in g/kg
    #
    #T_s=302
    #h=7
    #gamma_PH=6.2
    #zeta_T=16
    #p_s = 1000e2
    #
    #def findT(T_s, p):
    #    
    #    zeta = -h*(np.log(p/p_s))
    #    
    #    if (zeta < zeta_T):
    #        T = T_s - gamma_PH*zeta
    #    else:
    #        T = T_s - gamma_PH*zeta_T
    #    return T
    #    
    #qv_s = np.zeros(p.size)
    #    
    #for i, plev in enumerate(p):
    #    qv = wsat(findT(T_s, plev), p_s)
    #    qv_s[i] = qv*1e3    
    #    
    #qv_s = qv_s[interior]

    
    #if varname == 'W':
    #    axarr[1,].axvline(w_BL, color='r', alpha=0.5, label='simple model')
    #    axarr[0,].axvline(w_m, color='r', alpha=0.5, label='simple model')
    #if varname == 'QRAD':
    #    axarr[0,].axvline(W_m, color='r', label = 'simple model')
    #    axarr[1,].axvline(W_d, color='r', label = 'simple model')
    #    axarr[0,].set_xlim(-4, 3.5)
    #if varname == 'QV':
    #    axarr[1,].axvline(q_FA, color='r', alpha=0.5, label = 'simple model')
    #    axarr[0,].plot(qv_s, color='r', alpha=0.5, label = 'simple model')
    #
    

















