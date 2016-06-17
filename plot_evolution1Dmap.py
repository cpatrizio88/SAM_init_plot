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
fout = '/Users/cpatrizio/figures/SST302/SAM_aggr200days_12288km_64vert_ubarzero_1DMAPS/'

nc_in = glob.glob(fpath + '*4096x64*3000m*200days_302K.nc')[0]
nc_data= Dataset(nc_in)
varis = nc_data.variables

x = varis['x'][:]
t = varis['time'][:]

xx, tt = np.meshgrid(x, t)

#LWNT: OLR
#Prec: surface precip
varname = 'PW'
vari = varis[varname] 
tstep=int(12)
field = vari[:]
minvar= np.min(field)
maxvar= np.max(field)

plt.figure(1)
#plt.contour(xx/1e3, tt, field, np.linspace(minvar, maxvar, 12), colors='k', alpha=0.5)
#q = plt.quiver(xx[::8, ::8]/1e3, yy[::8, ::8]/1e3, usfc[i,::8,::8], vsfc[i,::8,::8], scale=500, alpha=0.8, zorder=1)
plt.contourf(xx/1e3, tt, field, np.linspace(minvar, maxvar, 24),cmap=cm.RdYlBu_r, zorder=0)
#p = plt.quiverkey(q, np.min(x)/1e3+30, np.max(y)/1e3+5, 15, "15 m/s",coordinates='data',color='k', alpha=0.8)
plt.xlabel('x (km)')
plt.ylabel('t (days)')
plt.title('{:s} [{:s}]'.format(varname, vari.units))
cb = plt.colorbar()
cb.set_label('({:s})'.format(vari.units))
plt.savefig(fout + '{:s}evolution_1D'.format(varname))
plt.close()

#test plotting 1-day average hovmoller plot.. 
#also test plotting at a single time (to compare with radial plots in 3D simulations)








