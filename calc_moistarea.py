from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
from thermolib.wsat import wsat
from thermolib.constants import constants
import SAM_init_plot.misc_fns 
from SAM_init_plot.misc_fns import fracclusterarea
from SAM_init_plot.block_fns import blockave2D

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (20, 10)})
matplotlib.rcParams.update({'lines.linewidth': 3})

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'

nc_in1 = glob.glob(fpath + '*256x256*3000m*day230to250*302K.nc')[0]
nc_in2 = glob.glob(fpath + '*512x512*3000m*day180to195*302K.nc')[0]
nc_in3 = glob.glob(fpath + '*1024x1024*3000m*day170to180*302K.nc')[0]


nc_in = [nc_in1, nc_in2, nc_in3]
domsize = [768, 1536, 3072]

#average overndays days
ndays=10
ntave2D = 24
nave = ndays*ntave2D

db=16


#threshold factor for calculating fractional area
#1.5 standard deviations above the mean PW
a=1


for i, nc_f in enumerate(nc_in):
    if nc_f == nc_in1:
        t=-1
    if nc_f == nc_in2:
        t=-1
    if nc_f == nc_in3:
        t=-1
    nc_data = Dataset(nc_f)
    varis2D = nc_data.variables
    x = varis2D['x'][:]
    y = varis2D['y'][:]
    times = varis2D['time'][:]
    daytimes=np.arange(0, times[-1]+1)
    #print 'domain size ({:3.0f} km)^2 fields averaged at day {:2.0f} over previous {:2.1f} days'.format(domsize[i], times[t], ndays)
    nx = x.size
    ny = y.size
    nxprime = nx / db
    nyprime = ny / db
    totpointsprime = nxprime*nyprime
    xx, yy = np.meshgrid(x, y)
    PW = varis2D['PW'][t-nave:t,:,:]
    W500 = varis2D['W500'][t-nave:t,:,:]
    P = varis2D['Prec'][t-nave:t,:,:]
    PWtave= np.mean(PW, axis=0)
    W500tave = np.mean(W500, axis=0)
    W500block = blockave2D(W500tave, db)
    Ptave = np.mean(P, axis=0)
    #fracA_Prec= fracclusterarea('Prec', varis2D, nave, t, a)
    fracA_PW = fracclusterarea('PW', varis2D, nave, t, a)
    W500crit = 0.01
   # fracA_W500  = fracclusterarea('W500', varis2D, nave, t, a, W500crit)
    totpoints=nx*ny
    fracA_W500 = len(W500block[W500block >= W500crit])/(1.*totpointsprime)
   

    fracA_Phat = fracclusterarea('PrecNEW', varis2D, nave, t)
    PWcrit = np.mean(PWtave) + a*np.std(PWtave)
    Pcrit = np.mean(Ptave) + a*np.std(Ptave)
 
    
    print 'moist edge PW threshold {:2.1f} (mm), mean PW = {:2.1f}, std PW = {:2.1f}'.format(PWcrit, np.mean(PWtave), np.std(PWtave))
    print 'moist edge precip threshold {:2.1f} (mm/day), mean precip = {:2.2f}, std precip = {:2.2f}, precip hat = {:2.2f}'.format(Pcrit, np.mean(Ptave), np.std(Ptave), np.mean(Ptave[Ptave > 0]))
    print 'moist edge W500 threshold {:3.3f} (m/s), mean W500 = {:3.3f}, std W500 = {:3.3f}'.format(W500crit, np.mean(W500tave), np.std(W500tave))
    #print 'domain size ({:d} km)^2, fractional area of convective region using precip threshold = {:3.1f}%'.format(domsize[i], fracA_Prec*100)
    print 'domain size ({:d} km)^2, fractional area of convective region using PW threshold = {:3.1f}%'.format(domsize[i], fracA_PW*100)
    print 'domain size ({:d} km)^2, fractional area of convective region using W500 threshold (W500 > {:3.4f}) = {:3.1f}%'.format(domsize[i], W500crit, fracA_W500*100)
    print 'domain size ({:d} km)^2, fractional area of convective region using precip bar/precip hat = {:3.1f}%'.format(domsize[i], fracA_Phat*100)
    r_m = np.sqrt(fracA_W500*domsize[i]**2)
    L = (1./np.sqrt(2))*domsize[i] #1/sqrt(2)*domsize is the max. radial distance possible from any point in a square domain
    print 'size of convective region using sigma w500 = {:4.3f} km'.format(r_m)
    print 'max radial distance from convective region center = {:4.3f} km'.format(L)
    print '---------------------------------------------------------------------'
    
    #plt.figure()
    #levels=np.arange(0, 90, 2)
    #cs = plt.contourf(xx/1e3, yy/1e3, PWtave, levels, cmap=cm.RdYlBu_r)
    ##plt.contour(xx/1e3, yy/1e3, W500tave, colors=['k'])
    #plt.xlabel('x (km)')
    #plt.ylabel('y (km)')
    #plt.title('PW (mm), {:d} km^2'.format(domsize[i]))
    #cb=plt.colorbar(cs)
    #cb.set_label('mm')
    #plt.show()
    #
    #plt.figure()
    #levels=np.arange(0, 55, 5)
    #cs = plt.contourf(xx/1e3, yy/1e3, Ptave, levels, cmap=cm.RdYlBu_r)
    #plt.contour(xx/1e3, yy/1e3, Ptave, levels=[Pcrit], colors=['k'])
    #plt.xlabel('x (km)')
    #plt.ylabel('y (km)')
    #plt.title('P (mm/day), {:d} km^2'.format(domsize[i]))
    #cb=plt.colorbar(cs)
    #cb.set_label('mm/day')
    #plt.show()

    
    
    



