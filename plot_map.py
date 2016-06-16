from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker
import matplotlib

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/figures/SST302/SAM_aggrday90to120_3072km_64vert_ubarzero_MAPS/'

nc_in = glob.glob(fpath + '*1024x1024*3000m*day90to120_302K.nc')[0]

nc_data= Dataset(nc_in)
varis = nc_data.variables
times = varis['time'][:]

x = varis['x'][:]
y = varis['y'][:]

xx, yy = np.meshgrid(x, y)
#LWNT: OLR
#Prec: surface precip
varname = 'LWNT'
vari = varis[varname]
tstep=int(12)
field = vari[:]
usfc = varis['USFC'][:]
vsfc = varis['VSFC'][:]
minvar= np.min(field)
maxvar= np.max(field)
#minvar = 20
#maxvar = 250
tstart=60*tstep
tstart=0

for i in np.arange(tstart, len(times), tstep):
    plt.figure()
    plt.contour(xx/1e3, yy/1e3, field[i,:,:], np.linspace(minvar, maxvar, 12), colors='k', alpha=0.5)
    #q = plt.quiver(xx[::8, ::8]/1e3, yy[::8, ::8]/1e3, usfc[i,::8,::8], vsfc[i,::8,::8], scale=500, alpha=0.8, zorder=1)
    plt.contourf(xx/1e3, yy/1e3, field[i,:,:], np.linspace(minvar, maxvar, 12),cmap=cm.RdYlBu_r, zorder=0)
    #p = plt.quiverkey(q, np.min(x)/1e3+30, np.max(y)/1e3+5, 15, "15 m/s",coordinates='data',color='k', alpha=0.8)
    plt.xlabel('x (km)')
    plt.ylabel('y (km)')
    plt.title('{:s} [{:s}] at t = {:3.2f} days'.format(varname, vari.units, times[i]))
    plt.colorbar()
    #cb = plt.colorbar(ticks=np.arange(10,90), 10)
    plt.savefig(fout + '{:s}map_day{:2.1f}.pdf'.format(varname,times[i]))
    plt.close()








