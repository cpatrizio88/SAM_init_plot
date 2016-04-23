import numpy as np
import site
site.addsitedir('/Users/cpatrizio/Dropbox/research/code/thermlib/')
site.addsitedir('/Users/cpatrizio/Dropbox/research/')
from theta import theta
from findT import *

p_s = 1013e2
T_s = 302

p_BL = 950e2
p_trop = 150e2

fpath =  '/Users/cpatrizio/SAM6.10.8/RCE_IDEAL/'
fname_in = 'lsf_TOGA'

lsf_in = np.loadtxt(fpath + fname_in, skiprows=0)

plevs = lsf_in[:,1]*1e2

levels=plevs.shape
lsf_out = np.ones((plevs.size, 7))
zlevs = -999*np.ones(plevs.shape)
T = np.zeros(plevs.shape)
thetat=np.zeros(plevs.shape)
qt = np.zeros(plevs.shape)
u_bar = np.zeros(plevs.shape)
v_bar = np.zeros(plevs.shape)
w_bar = np.zeros(plevs.shape)

lsf_out[:,0] = zlevs
lsf_out[:,1] = plevs/100.
lsf_out[:,2] = thetat
lsf_out[:,3] = qt
lsf_out[:,4] = u_bar
lsf_out[:,5] = v_bar
lsf_out[:,6] = w_bar

day0=0

fname_out = 'lsf_ideal'
head= 'z[m] p[mb] tpls[K/s] qls[kg/kg/s] uls_hor vls_hor wls[m/s]'

np.savetxt(fpath + fname_out, lsf_out, delimiter=' ', header=head, fmt='%3.4f')




    
    