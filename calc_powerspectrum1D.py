from netCDF4 import Dataset
import numpy as np
import matplotlib
import glob
import matplotlib.pyplot as plt


matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/figures/SST302/SAM_aggr200days_12288km_64vert_ubarzero_1DFFT/'

nc_in = glob.glob(fpath + '*4096x64*3000m*200days_302K.nc')[0]
nc_data= Dataset(nc_in)
varis = nc_data.variables

Lx = 12288
x=varis['x'][:]
t=varis['time'][:]

Ls_kave = np.zeros(t.size)

n = len(x)

#####USER EDIT: change to look at different variables and/or time
varname='PW'

field = varis[varname][:]

#calculate for different times
for i in np.arange(t.size):
    field_t = field[i, :]
    #A contains fourier coefficients
    A = np.fft.fft(field_t, n)
    #phi is the power spectrum
    phi = np.abs(A)**2
    #Nyquist wavenumber
    k_nyq = n/2
    
    phipos = phi[1:k_nyq]
    
    k = np.arange(1, k_nyq)
    #average wavenumber
    kave = np.sum(np.multiply(k, phipos))/np.sum(phipos)
    #average length 
    Ls_kave[i] = Lx/kave


plt.figure(1)
plt.plot(t, Ls_kave)
plt.xlabel('time (days)')
plt.ylabel(r'$L_{{k}}$ (km)')
plt.title(r'power-weighted average wavenumber of spatial variability of {0}, $L_k$, domain size {1} km'.format(varname, Lx))
plt.savefig(fout + '{0}corrlength.pdf'.format(varname))
plt.close()





