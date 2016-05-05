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
fout = '/Users/cpatrizio/figures/SAM RCE 100 days figs/'

nc_in = glob.glob(fpath + '*100days.nc')[0]

nc_data= Dataset(nc_in)
vari = nc_data.variables
times = vari['time'][:]

x = vari['x'][:]
y = vari['y'][:]

xx, yy = np.meshgrid(x, y)

varname = 'PW'
tstep=int(24*8)
field = vari[varname][:]

for t in times[::tstep]:
    plt.figure()
    plt.contourf(xx/1e3, yy/1e3, field[t,:,:], 30, cmap=cm.RdYlBu)
    plt.xlabel('x (km)')
    plt.ylabel('y (km)')
    plt.title('{:s} at t = {:3.2f}'.format(varname, t))
    plt.colorbar()
    plt.show()












