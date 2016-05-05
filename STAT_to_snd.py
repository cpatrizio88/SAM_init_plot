from netCDF4 import Dataset
import glob
import numpy as np

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
snd_path =  '/Users/cpatrizio/SAM6.10.8/RCE_IDEAL/'

snd_inname = 'snd_ideal'

nc_in = glob.glob(fpath + '*100days.nc')[0]
nc_data = Dataset(nc_in)
nc_vars = nc_data.variables

z = nc_vars['z'][:]
p = nc_vars['p'][:]

qv_STAT = nc_vars['QV'][:]
theta_STAT = nc_vars['THETA'][:]

qv_RCE = qv_STAT[-1,:]
theta_RCE = theta_STAT[-1,:]


snd_in = np.loadtxt(snd_path + snd_inname, skiprows=0)

snd_out = np.ones((p.size, 6))
zlevs = -999*np.ones(p.shape)
q = np.zeros(p.shape)
T = np.zeros(p.shape)
thetas=np.zeros(p.shape)
u = np.zeros(p.shape)
v = np.zeros(p.shape)

snd_out[:,0] = z
snd_out[:,1] = p
snd_out[:,2] = theta_RCE
snd_out[:,3] = qv_RCE
snd_out[:,4] = u
snd_out[:,5] = v

fname_out = 'snd_RCE'
head='z[m] p[mb] T[K] q[g/kg] u[m/s] v[m/s]'

np.savetxt(snd_path + fname_out, snd_out, delimiter=' ', header=head, fmt='%3.4f')