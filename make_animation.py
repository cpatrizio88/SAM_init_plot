from netCDF4 import Dataset
import matplotlib.animation as animation
import types
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
fout = '/Users/cpatrizio/figures/SST302/1536km_SAM_aggrday90to130_64vert_ubarzero_MAPS/'

nc_in = glob.glob(fpath + '*512x512*3000m*day90to130*302K.nc')[0]

nc_data= Dataset(nc_in)
varis = nc_data.variables
times = varis['time'][:]

x = varis['x'][:]
y = varis['y'][:]

xx, yy = np.meshgrid(x, y)
#LWNT: OLR
#Prec: surface precip
varname = 'PW'
vari = varis[varname]
tstep=1
field = vari[:]
usfc = varis['USFC'][:]
vsfc = varis['VSFC'][:]
minvar= np.min(field)
maxvar= np.max(field)
#minvar = 20
#maxvar = 250
tstart=60*tstep
tstart=0

fig = plt.figure()
plt.xlabel('x (km)')
plt.ylabel('y (km)')

plt.contourf(xx/1e3, yy/1e3, field[0,:,:], np.linspace(minvar, maxvar, 12),cmap=cm.RdYlBu_r, zorder=0)
plt.colorbar()

def setvisible(self,vis):
    for c in self.collections: c.set_visible(vis)
def setanimated(self,ani):
    for c in self.collections: c.set_animated(ani)

ims = []

for i in np.arange(tstart, len(times), tstep):
    print(i)
    ax = plt.gcf().gca()
    title = ax.set_title('{:s} [{:s}] at t = {:3.2f} days'.format(varname, vari.units, times[i]))
    im2 = plt.contourf(xx/1e3, yy/1e3, field[i,:,:], np.linspace(minvar, maxvar, 12),cmap=cm.RdYlBu_r, zorder=0)
    im1 = plt.contour(xx/1e3, yy/1e3, field[i,:,:], np.linspace(minvar, maxvar, 12), colors='k', alpha=0.5)
    im1.set_visible = types.MethodType(setvisible,im1)
    im1.set_animated = types.MethodType(setanimated,im1)
    im1.axes = plt.gca()
    im1.figure=fig
    im2.set_visible = types.MethodType(setvisible,im2)
    im2.set_animated = types.MethodType(setanimated,im2)
    im2.axes = plt.gca()
    im2.figure=fig
 
    ims.append([im1, im2, title])
    #ims.append([im2])
    
ani = animation.ArtistAnimation(fig, ims, interval=10, blit=False, repeat_delay=1000)
FFwriter = animation.FFMpegWriter()
fout_ani =  '/Users/cpatrizio/animations/'
print 'saving animation'
ani.save(fout_ani + '1536kmPWanimation_test.mp4', writer = FFwriter)





    
    









