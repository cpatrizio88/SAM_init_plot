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
matplotlib.rcParams.update({'figure.figsize': (20, 10)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 16})

plt.style.use('seaborn-white')

fpath2D =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpathSTAT = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'

foutdata = '/Users/cpatrizio/data/SST302/'

fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr130days_64vert_ubarzero_MOISTDRYPROFS/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to140_64vert_ubarzero_MOISTDRYPROFS/'

nc_inSTAT = glob.glob(fpathSTAT + '*256x256*3000m*130days*302K.nc')[0]
nc_in2D = glob.glob(fpath2D + '*256x256*3000m*130days*302K.nc')[0]
nc_in3D = glob.glob(fpath3D + '*256x256*3000m*130days*302K.nc')[0]

#nc_inSTAT = glob.glob(fpathSTAT + '*1024x1024*3000m*150days*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*1024x1024*3000m*day140to150*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*1024x1024*3000m*day140to150*302K.nc')[0]

domsize=768

nave=5
ntave2D=24
ntave3D=4
t2 = -35*ntave2D
t3 = -35*ntave3D

aveperiod2D = nave*ntave2D
aveperiod3D = nave*ntave3D

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

#varnames = ['QRAD', 'W', 'U', 'QN', 'QV']
#varnames=['QV', 'TABS']
varnames=['QV', 'QN']
#varnames=['TABS']

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


#EDIT HIS TO CHANGE VARIABLE FOR MOIST REGION THRESHOLD
mvarname = 'PW'
mvar = varis2D[mvarname]
mfield = mvar[t2-aveperiod2D:t2,:,:]
mfield_t= np.mean(mfield, axis=0)

W500 = varis2D['W500'][t2-aveperiod2D:t2,:,:]
W500_t = np.mean(W500, axis=0)
W500crit=0.1

#EDIT THIS TO CHANGE HOW MANY STD ABOVE MEAN FOR MOIST REGION THRESHOLD
a=1.5

mfieldcrit = np.mean(mfield) + a*np.std(mfield)
moist_points = mfield_t >= mfieldcrit
dry_points = mfield_t <= mfieldcrit 

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

U_rfname = glob.glob(foutdata + '*{:d}*'.format(domsize))[0]

U_r = np.load(foutdata + U_rfname)


for varname in varnames:
    if varname == 'URADIAL':
        field_tave = U_r
        units = 'm/s'
    else: 
        print 'loading', varname
        vari = varis3D[varname]
        units = vari.units.strip()
        field = varis3D[varname][t3-aveperiod3D:t3,:,:,:]
        field_tave = np.mean(field, axis=0)
        #calculate relative humidity 
        #if varname == 'TABS':
        #    varname = 'RH'
        #    QV = varis3D['QV'][t3-aveperiod3D:t3,:,:,:]
        #    QV_tave = np.mean(QV, axis=0)
        #    RH_tave = np.zeros(QV_tave.shape)
        #    for i, plev in enumerate(p):
        #        wsat_tave = 1000*wsat(field_tave[i,:,:], plev) #convert to g/kg
        #        RH_tave[i,:,:] = 100*(QV_tave[i,:,:]/wsat_tave) #convert to percent
        #    field_tave = RH_tave
        
        moist_field = field_tave[:,moist_points]
        dry_field = field_tave[:,dry_points]
        
        moist_prof = np.mean(moist_field, axis=1)
        dry_prof = np.mean(dry_field, axis=1)
        
        if varname == 'RH':
            units='%'
        
        #simple model results   
        w_BL = -.0066
        w_m = 0.0334
        
        W_m = 3
        W_d = -4.5
        
        q_FA = 0.000175*1e3 #in g/kg
        
        T_s=302
        h=7
        gamma_PH=6.2
        zeta_T=16
        p_s = 1000e2
        
        def findT(T_s, p):
            
            zeta = -h*(np.log(p/p_s))
            
            if (zeta < zeta_T):
                T = T_s - gamma_PH*zeta
            else:
                T = T_s - gamma_PH*zeta_T
            return T
            
        qv_s = np.zeros(p.size)
            
        for i, plev in enumerate(p):
            qv = wsat(findT(T_s, plev), p_s)
            qv_s[i] = qv*1e3
            
        interior = np.bitwise_and(z > z_BL, z < z_t)
        zint = z[interior]
        moist_prof = moist_prof[interior] 
        dry_prof = dry_prof[interior]
        
        moist_profave = np.mean(moist_prof)
        dry_profave = np.mean(dry_prof)
        
        qv_s = qv_s[interior]
        
        plt.figure()
        f, axarr = plt.subplots(2,1)
        plt.suptitle(r'{:s} ({:s}), day {:3.0f} to {:3.0f} average, domain size = {:d} km$^2$, convective region {:s} threshold = {:2.1f} {:s} '.format(varname, units, t3D[0], t3D[-1], domsize, mvarname, mfieldcrit, mvar.units.strip()))
        axarr[0,].plot(moist_prof, zint/1e3, 'k')
        axarr[0,].set_title('{:s}, convection region'.format(varname))
        if varname == 'W':
            axarr[1,].axvline(w_BL, color='r', alpha=0.5, label='simple model')
            axarr[0,].axvline(w_m, color='r', alpha=0.5, label='simple model')
        if varname == 'QRAD':
            axarr[0,].axvline(W_m, color='r', label = 'simple model')
            axarr[1,].axvline(W_d, color='r', label = 'simple model')
            axarr[0,].set_xlim(-4, 3.5)
        if varname == 'QV':
            axarr[1,].axvline(q_FA, color='r', alpha=0.5, label = 'simple model')
            axarr[0,].plot(qv_s, color='r', alpha=0.5, label = 'simple model')
        axarr[0,].set_xlabel('{:s} ({:s})'.format(varname, units))
        axarr[0,].set_ylabel('z (km)')
        axarr[0,].set_ylim(z_BL/1e3, z_t/1e3)
        if (varname != 'QV'):
            axarr[0,].axvline(moist_profave, color='k', alpha=0.5, label='vertical average')
        axarr[1,].plot(dry_prof, zint/1e3, 'k')
        axarr[1,].set_title('{:s}, dry region'.format(varname))
        axarr[1,].set_xlabel('{:s} ({:s})'.format(varname, units))
        axarr[1,].set_ylabel('z (km)')
        axarr[1,].set_ylim(z_BL/1e3, z_t/1e3)
        axarr[1,].axvline(dry_profave, color='k', alpha=0.5, label='vertical average')
    
        plt.legend()
        plt.xlabel('{:s} ({:s})'.format(varname, vari.units.strip()))
        plt.savefig(fout + '{:s}dryconvprof_day{:3.0f}to{:3.0f}.pdf'.format(varname, t3D[0], t3D[-1]))
        plt.close()
        
        varBL_moist = moist_prof[BLi]
        varBL_dry = dry_prof[BLi]
        
        print '{:s} at convective region boundary layer top, z = {:3.1f}, is: {:4.4f} ({:s})'.format(varname, z_BL, varBL_moist, units)
        print '{:s} at dry region boundary layer top, z = {:3.1f}, is: {:4.4f} ({:s})'.format(varname, z_BL, varBL_dry, units)
        
    
    
    
    
















