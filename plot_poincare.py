import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import matplotlib
import matplotlib.cm as cm

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.figsize': (20, 13)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 16})

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
fout = '/Users/cpatrizio/Google Drive/figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STAT/'

nc_in1 = glob.glob(fpath + '*256x256*3000m*180days_302K.nc')[0]
#nc_in2 = glob.glob(fpath + '*512x512*3000m*140days_302K.nc')[0]
#nc_in3 = glob.glob(fpath + '*1024x1024*3000m*150days_302K.nc')[0]

nc_data = Dataset(nc_in1)

nc_vars = nc_data.variables

nave2D=24

tinit = 60*nave2D
tend = -1

tperiod = 60*nave2D

t = nc_vars['time'][tinit:tend]


PW = nc_vars['PW'][tinit:tend]

PWprev = PW[:-1]
PWnext = PW[1:]

varname = 'LWNT'
var = nc_vars[varname]
tseries = nc_vars[varname][tinit:tend]

n=24

tseriessmooth = moving_average(tseries,n)

PWsmooth = moving_average(PW, n)

PWprevsmooth = moving_average(PWprev, n)
PWnextsmooth = moving_average(PWnext, n)

#poincare = zip(PW, LWNT)

plt.figure()
#ax = plt.gcf().gca()
fig, axarr = plt.subplots(2,1)
axarr[0,].plot(PWprev[:tperiod], PWnext[:tperiod], 'yx', label='day {:3.0f} to {:3.0f}'.format(t[0], t[tperiod]))
#plt.plot(PWprev[tperiod:], PWnext[tperiod:], 'x')
#plt.plot(PWprev[0], PWnext[0], 'yo', markersize=10)
#plt.plot(PWprev[tperiod], PWnext[tperiod], 'bo', markersize=10)
#plt.plot(PWprev[-1], PWnext[-1], 'ro', markersize=10)
axarr[1,].plot(PWprev[tperiod:], PWnext[tperiod:], 'bx', label='day {:3.0f} to {:3.0f}'.format(t[tperiod], t[-1]))
plt.xlabel('PW (mm)')
plt.ylabel(r'PW (mm)')
plt.suptitle('Poincare plot, day {:3.0f} to day {:3.0f}'.format(t[0], t[-1]))
plt.legend()
plt.show()

plt.figure()
#ax = plt.gcf().gca()
fig, axarr = plt.subplots(2,1)
axarr[0,].plot(PWprevsmooth[:tperiod], PWnextsmooth[:tperiod], 'yx', label='day {:3.0f} to {:3.0f}'.format(t[0], t[tperiod]))
axarr[1,].plot(PWprevsmooth[tperiod:], PWnextsmooth[tperiod:], 'bx', label='day {:3.0f} to {:3.0f}'.format(t[tperiod], t[-1]))
#plt.plot(PWprevsmooth[0], PWnextsmooth[0], 'yo', markersize=10, label='day {:3.0f}'.format(t[0]))
#plt.plot(PWprevsmooth[tperiod], PWnextsmooth[tperiod], 'bo', markersize=10, label='day {:3.0f}'.format(t[tperiod]))
#plt.plot(PWprevsmooth[-1], PWnextsmooth[-1], 'ro', markersize=10, label='day {:3.0f}'.format(t[-1]))
plt.xlabel('PW (mm)')
plt.ylabel(r'PW (mm)')
plt.suptitle('Poincare plot, smoothed n = {:d} hour window, day {:3.0f} to day {:3.0f}'.format(n, t[0], t[-1]))
plt.legend()
plt.show()

plt.figure()
plt.plot(PWsmooth[:tperiod], tseriessmooth[:tperiod], 'x', label='day {:3.0f} to {:3.0f}'.format(t[0], t[tperiod]))
plt.plot(PWsmooth[tperiod:], tseriessmooth[tperiod:], 'x', label='day {:3.0f} to {:3.0f}'.format(t[tperiod], t[-1]))
plt.plot(PWsmooth[0], tseriessmooth[0], 'yo', markersize=10, label='day {:3.0f}'.format(t[0]))
plt.plot(PWsmooth[tperiod], tseriessmooth[tperiod], 'bo', markersize=10, label='day {:3.0f}'.format(t[tperiod]))
plt.plot(PWsmooth[-1], tseriessmooth[-1], 'ro', markersize=10, label='day {:3.0f}'.format(t[-1]))
plt.title('Poincare plot, smoothed n = {:d} hour window, day {:3.0f} to day {:3.0f}'.format(n, t[0], t[-1]))
plt.xlabel('PW (mm)')
plt.ylabel('{:s} ({:s})'.format(varname, var.units.strip()))
plt.legend()
plt.show()

plt.figure()
plt.plot(PWprevsmooth[:tperiod], label='day {:3.0f} to {:3.0f}'.format(t[0], t[tperiod]))
plt.plot(PWprevsmooth[tperiod:], label='day {:3.0f} to {:3.0f}'.format(t[tperiod], t[-1]))
plt.title('PW smoothed with n = {:d} hour window'.format(n))
plt.xlabel('time (days)')
plt.ylabel('PW (mm)')
plt.legend()
plt.show()

plt.figure()
plt.plot(tseriessmooth[:tperiod], label='day {:3.0f} to {:3.0f}'.format(t[0], t[tperiod]))
plt.plot(tseriessmooth[tperiod:], label='day {:3.0f} to {:3.0f}'.format(t[tperiod], t[-1]))
plt.title('{:s} smoothed with n = {:d} hour window'.format(varname, n))
plt.xlabel('time (days)')
plt.ylabel('{:s} ({:s})'.format(varname, var.units.strip()))
plt.legend()
plt.show()









