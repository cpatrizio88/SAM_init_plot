from netCDF4 import Dataset
import numpy as np
import matplotlib
import glob
import matplotlib.pyplot as plt


matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})

plt.style.use('seaborn-white')

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_2D/'
fout = '/Users/cpatrizio/figures/SST302/SAM_aggrday90to125_1536km_64vert_ubarzero_2DFFT/'

nc_in = glob.glob(fpath + '*512x512*3000m*day90to125*302K.nc')[0]
nc_data= Dataset(nc_in)
varis = nc_data.variables

Lx = 768
x=varis['x'][:]
y=varis['y'][:]
nx=len(x)
ny = len(y)

k_nyq = nx/2
l_nyq = ny/2
k = np.arange(1, k_nyq)
l = np.arange(1, l_nyq)
#convert to 2D for later calculations
l = np.tile(l, len(k)).reshape(len(k), len(l))

times = np.arange(0, 24*120)

Ls_kave = np.zeros(times.shape)

#####USER EDIT: change to look at different variables and/or time
varname='PW'

field = varis[varname][:]

#calculate for different times
for i,t in enumerate(times):
    field_t = field[t, :]
    #A contains fourier coefficients
    A = np.fft.fft2(field_t)
    #phi is the power spectrum
    phi = np.power(np.abs(A), 2)
    #Nyquist wavenumber
    phipos = phi[1:k_nyq, 1:l_nyq]
    nume=np.sum(np.multiply(k, np.sum(np.multiply(l, phipos), axis=1)), axis=0)
    #average wavenumber
    kave = nume/np.sum(phipos)
    #average length 
    Ls_kave[i] = Lx/kave

plt.figure(1)
plt.plot(times/24., Ls_kave)
plt.xlabel('time (days)')
plt.ylabel(r'correlation length, $L_{{k}}$ (km)')
plt.title(r'correlation length of {0}, domain size {1} km$^2$'.format(varname, Lx))
plt.savefig(fout + '{0}corrlength.pdf'.format(varname))
plt.close()





