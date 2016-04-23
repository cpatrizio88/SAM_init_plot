import numpy as np
import site
site.addsitedir('/Users/cpatrizio/Dropbox/research/code/thermlib/')
site.addsitedir('/Users/cpatrizio/Dropbox/research/')
from findT import *
from theta import *
import matplotlib.pyplot as plt

fpath =  '/Users/cpatrizio/SAM6.10.8/RCE_IDEAL/'
fname = 'snd_TOGA'

snd = np.loadtxt(fpath + fname, skiprows=0)
p_s = 1013e2
p_BL = 950e2
p_trop = 150e2
gamma = 6.5

p = snd[:,1]
theta_TOGA = snd[:,2]
qv = snd[:,3]

T_constgamma = np.zeros(p.shape)
theta_constgamma = np.zeros(p.shape)
T_s = theta_TOGA[0]

    
for i, pval in enumerate(p):
    T_constgamma[i] = findT(T_s, pval*100, p_s, p_trop, gamma)
    theta_constgamma[i] = theta(T_constgamma[i], pval*100)
    
plt.figure(1)
plt.plot(theta_constgamma, p, 'kx', label='gamma = {0} K km^-1'.format(gamma))
plt.plot(theta_TOGA, p,'rx', label='theta_TOGA')
plt.xlabel('potential temperature (K)')
plt.ylabel('pressure (hPa)')
plt.gca().invert_yaxis()
plt.legend()
plt.show()

q_BL = 0.017
q_FA = 0.001
q_strat = 1e-6
qv_ideal = np.zeros(qv.shape)

qv_ideal[p*1e2 >= p_BL] = q_BL
qv_ideal[p*1e2 < p_BL] = q_FA
qv_ideal[p*1e2 < p_trop] = q_strat 

qv_ideal = qv_ideal*1e3

plt.figure(2)
plt.plot(qv, p, 'rx', label='qv_TOGA')
plt.plot(qv_ideal, p, 'kx', label='qv_ideal')
plt.xlabel('mixing ratio (g/kg)')
plt.ylabel('pressure (hPa)')
plt.gca().invert_yaxis()
plt.legend()
plt.show()

    
    
    



