from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
import matplotlib.cm as cm
import SAM_init_plot.block_fns
import SAM_init_plot.misc_fns
import thermolib
from SAM_init_plot.misc_fns import radprof3D
from SAM_init_plot.block_fns import blockave2D, blockave3D, blockxysort2D
from thermolib.constants import constants
from thermolib.wsat import wsat

c = constants()

matplotlib.rcParams.update({'font.size': 24})
matplotlib.rcParams.update({'figure.figsize': (22, 14)})
matplotlib.rcParams.update({'legend.fontsize': 22})
plt.style.use('seaborn-white')

fpath3D =  '/Users/cpatrizio/SAM6.10.8/OUT_3D/'
fpath2D = '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
foutdata = '/Users/cpatrizio/data/SST302/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/768km_SAM_aggr250days_64vert_ubarzero_RADIALXSECTIONNEW/'
#fout = '/Users/cpatrizio/Google Drive/figures/SST302/1536km_SAM_aggrday90to195_64vert_ubarzero_RADIALXSECTIONNEW/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/3072km_SAM_aggrday110to190_64vert_ubarzero_RADIALXSECTIONNEW/'

#nc_in3D = glob.glob(fpath3D + '*256x256*3000m*day230to250*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*256x256*3000m*day230to250*302K.nc')[0]
#nc_in3D = glob.glob(fpath3D + '*512x512*3000m*day180to195*302K.nc')[0]
#nc_in2D = glob.glob(fpath2D + '*512x512*3000m*day180to195*302K.nc')[0]
nc_in3D = glob.glob(fpath3D + '*1024x1024*3000m*day170to180*302K.nc')[0]
nc_in2D = glob.glob(fpath2D + '*1024x1024*3000m*day170to180*302K.nc')[0]


#domsize=768
#domsize=1536
domsize=3072

nave=3
ntave2D=24
ntave3D=4
t2 = -35*ntave2D
t3 = -35*ntave3D

t2=-1
t3=-1

aveperiod3D = nave*ntave3D
aveperiod2D = nave*ntave2D

print 'loading .nc files'
nc_data3D = Dataset(nc_in3D)
nc_data2D = Dataset(nc_in2D)
varis3D = nc_data3D.variables
varis2D = nc_data2D.variables

print 'loading 2D variables'
nc_data2D = Dataset(nc_in2D)
varis2D = nc_data2D.variables

print 'loading 3D variables'
nc_data3D = Dataset(nc_in3D)
varis3D = nc_data3D.variables

PW = varis2D['PW'][t2-aveperiod2D:t2,:,:]
PW_tave = np.mean(PW, axis=0)

U = varis3D['U'][t3-aveperiod3D:t3,:,:,:]
U_tave = np.mean(U, axis=0)
V = varis3D['V'][t3:aveperiod3D:t3,:,:,:]
V_tave = np.mean(V, axis=0)

t2D = varis2D['time'][t2-aveperiod2D:t2]
t3D = varis3D['time'][t3-aveperiod3D:t3]
x = varis3D['x'][:]
y = varis3D['y'][:]
z = varis3D['z'][:]
p = varis3D['p'][:]
p = p*1e2

xx, yy = np.meshgrid(x, y)
times = np.arange(t3D[0], np.max(t3D))
db=16

print 'calulating max PW, blocked'
#PW_t = np.mean(PW_tave[t-nave:t,:,:], axis=0)
#PW_t = np.mean(PW[t-nave2D:t,:,:], axis=0)
PW_tave = np.mean(PW, axis=0)

PW_blocked = blockave2D(PW_tave, db)
PWxy_sorted = blockxysort2D(PW_tave, xx, yy, db)

PWsort = PWxy_sorted.keys()
PWxycoords = PWxy_sorted.values()

mcenter = PWxycoords[-1]

#calculate distance from convective center, taking into account periodic boundaries
width = xx[0,-1]
delx = np.abs(xx - mcenter[0])
dely = np.abs(yy - mcenter[1])
delx[delx > width/2] = width - delx[delx > width/2]
dely[dely > width/2] = width - dely[dely > width/2]
r = np.sqrt(np.power(delx, 2) + np.power(dely, 2))

gradry, gradrx = np.gradient(r)

delx = delx*np.sign(gradrx)
dely = dely*np.sign(gradry)

delx = delx/r
delx3D = np.ones(U_tave.shape)
delx3D[:,:,:] = delx 
#delx3D = delx3D.T
delx3D = delx3D.flatten()

#delx = delx.flatten()
dely = dely/r
dely3D = np.ones(U_tave.shape)
#dely = dely.flatten()
#dely3D = dely3D.T
dely3D = dely3D.flatten()



rhat = zip(delx3D, dely3D)

Uflat = U_tave.flatten()
Vflat = V_tave.flatten()

Uvec = zip(Uflat, Vflat)

U_r = np.einsum('ij,ij->i', rhat, Uvec)


#u_r = U_r[:,0]
#v_r = U_r[:,1]

#U_r = np.sqrt(u_r**2 + v_r**2)
U_r = U_r.reshape(z.size, x.size, y.size)

np.save(foutdata + '{:d}km_URADIAL_day{:3.0f}to{:3.0f}'.format(domsize, t3D[0], t3D[-1]), U_r)

#save radial profile somewhere?


#plot radial profile

binwidth=5e3
rbins = np.arange(0, domsize*1e3, binwidth)

zedges = np.zeros(z.shape)
zedges[1:] = (z[1:] + z[:-1])/2.

nbins=[rbins, zedges]

redges, zedges, fieldmeans = radprof3D(U_r, xx, yy, z, mcenter, nbins=nbins)

rbin_centers = (redges[1:] + redges[:-1])/2.
zbin_centers = (zedges[1:] + zedges[:-1])/2.

#rr has shape (zedges.size, redges.size) = nbins.T
rr, zz = np.meshgrid(rbin_centers, zbin_centers)

fieldmeans = np.transpose(fieldmeans)

vmin=np.ma.masked_invalid(fieldmeans).min()
vmax=np.ma.masked_invalid(fieldmeans).max()
#EDIT: COLORBAR FORMATTING. IF NOT MONOTONIC, SWITCH TO DIVERGING COLORBAR, OTHERWISE
#USE SEQUENTIAL COLORBAR. SET COLOR MAPPING LIMITS MANUALLY. 
if (vmin < 0 and vmax > 0):
    if np.abs(vmin) < np.abs(vmax):
        vmin = vmin - (vmax - np.abs(vmin))
    else:
        vmax = vmax + (np.abs(vmin) - vmax)
    cmap = cm.RdBu_r
else:
    cmap = cm.YlGnBu
fig = plt.figure()
units = 'm/s'
varname = 'URADIAL'
titlename = r'$u_r$'
#vmin=3
#vmax=-3
#cmap = cm.RdYlBu_r
ax = fig.gca()
plt.pcolormesh(rr/(domsize*1e3), zz/1e3, fieldmeans, vmin=vmin, vmax=vmax, cmap=cmap)
#plt.xlabel('radial distance (km)')
plt.xlabel(r'$\hat{r}$')
plt.ylabel('z (km)')
cb=plt.colorbar()
cb.set_label(units)
ax.set_ylim(0,26)
ax.set_xlim([0, (1/np.sqrt(2))])
#ax.set_xlim([0, (1/np.sqrt(2))*domsize])
plt.title('{:s}, day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$'.format(titlename, t3D[0], t3D[-1], domsize))
#plt.title('{:s} ({:s}), day {:3.0f} to {:3.0f} average, domain size = ({:d} km)$^2$, x bin width = {:2.2f}'.format(varname, vari.units.strip(), t3D[t-nave*ntave3D], t3D[t], domsize, 1./nbins[0]))
plt.savefig(fout + '{:s}radialxsectiondist_day{:3.0f}to{:3.0f}.pdf'.format(varname, t3D[0], t3D[-1]))
#plt.savefig(fout + '{:s}radialxsection_day{:3.0f}to{:3.0f}.pdf'.format(varname, t3D[t-nave*ntave3D], t3D[t]))
plt.close()

plt.figure(2)

skip = (slice(None, None, 5), slice(None, None, 5))

rxhat = delx.reshape(x.size, y.size)
ryhat = dely.reshape(x.size, y.size)
#
##fig, ax = plt.subplots()
##im = ax.imshow(r/1e3, extent=[x.min()/1e3, x.max()/1e3, y.min()/1e3, y.max()/1e3], cmap = cm.viridis)
##ax.quiver(xx[skip]/1e3, yy[skip]/1e3, rxhat[skip]/1e3, ryhat[skip]/1e3)
##
###fig.colorbar(im)
###ax.set(aspect=1, title='r hat')
##plt.show()

#test = np.random.rand(x.size, y.size)
#
#plt.figure(3)
#plt.contourf(xx/1e3, yy/1e3, rxhat, 10, cmap = cm.inferno)
##plt.title('delx')
#plt.show()
#
#
#plt.figure(4)
#plt.contourf(xx/1e3, yy/1e3, ryhat, 10, cmap = cm.inferno)
##plt.title('dely')
#plt.show()
#
##plt.figure(5)
##plt.contourf(xx/1e3, yy/1e3, test, 5, linewidth=5)
##plt.contourf(xx/1e3, yy/1e3, test, 10, cmap=cm.jet)
##plt.show()
#
#plt.figure(6)
#plt.contourf(xx/1e3, yy/1e3, U_r[20,:,:], cmap=cm.RdYlBu)
#plt.show()
#
#plt.figure(7)
#plt.contourf(xx/1e3, yy/1e3, U_tave[20,:,:], cmap=cm.RdYlBu)
#plt.title('U')
#plt.show()
#
#plt.figure(7)
#plt.contourf(xx/1e3, yy/1e3, V_tave[20,:,:], cmap=cm.RdYlBu)
#plt.title('V')
#plt.show()









#test einsum. the follow works for dot product of a list of 2D vectors
#r = np.array([[1,2], [2,3], [3,4]])
#v = np.array([[4,5], [6,7], [8,9]])
#v_r = np.einsum('ij,ij->ij', r,v)




