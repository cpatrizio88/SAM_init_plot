import numpy as np
import site
import matplotlib.pyplot as plt
site.addsitedir('/Users/cpatrizio/Dropbox/research/code/thermlib/')
site.addsitedir('/Users/cpatrizio/Dropbox/research/')
from findT import *
from wsat import *
import matplotlib


matplotlib.rcParams.update({'font.size': 28})
matplotlib.rcParams.update({'figure.figsize': (18, 10)})
matplotlib.rcParams.update({'lines.linewidth': 1})
matplotlib.rcParams.update({'legend.fontsize': 22})

p_0=1000e2 #reference pressure (Pa)
Rd = 287.04   # Gas constant of dry air [J kg^-1 K^-1]
Cpd = 1005.7  # Heat capacity at constant pressure for dry air [J kg^-1 K^-1]

#   potential temperature [K] from pressure P [Pa] and temperature T [K]
def theta(p, T): return T*(p_0/p)**(Rd/Cpd)

def Tinvert(p, theta): return theta*(p/p_0)**(Rd/Cpd)

p_s = 1015e2


SSTs = [302]

for T_s in SSTs:
    p_BL = 950e2
    p_trop = 100e2
    
    fpath =  '/Users/cpatrizio/SAM6.10.8/RCE_IDEAL/'
    fname_in = 'snd_TOGA'
    
    
    snd_in = np.loadtxt(fpath + fname_in, skiprows=0)
    
    plevs = snd_in[:,1]*1e2
    
    q_toga = snd_in[:,3] #water vapor profile from TOGA (g/kg)
    theta_toga = snd_in[:,2] #theta profile from TOGA (K)
    RH_toga = np.zeros(plevs.shape)
    T_toga = np.zeros(plevs.shape)
    RH_new = np.zeros(plevs.shape)
    
    #get RH humidity profile from TOGA
    for i, p in enumerate(plevs):
        T_toga[i] = Tinvert(p, theta_toga[i])
        RH_toga[i] = q_toga[i]/(wsat(T_toga[i], p)*1e3)
        if (p <= p_BL) and (p > p_trop):
            RH_new[i] = 0.1 + 0.65*np.exp(6*(p - p_BL)/p_trop)
        elif p > p_BL:
            RH_new[i] = 0.75
        else:
            RH_new[i] = 0.01
    
        
    snd_out = np.ones((plevs.size, 6))
    zlevs = -999*np.ones(plevs.shape)
    q = np.zeros(plevs.shape)
    q_humid = np.zeros(plevs.shape)
    T = np.zeros(plevs.shape)
    thetas=np.zeros(plevs.shape)
    u = np.zeros(plevs.shape)
    v = np.zeros(plevs.shape)
    
    gamma_PH = 6.5 #lapse rate (K km^-1)
    T_pierrehumb = np.zeros(plevs.shape)
    
    #idealized temperature profile from Pierrehumbert (1995)
    for i, p in enumerate(plevs):
        T[i] = findT(T_s, p, p_s, p_trop, gamma_PH) 
        #set moisture profile to same RH as toga profile
        q[i] = RH_new[i]*wsat(T[i], p)
        q_humid[i] = RH_toga[i]*wsat(T[i], p)
        thetas[i] = theta(p, T[i])
    
    #convert to g/kg  
    q = q*1e3
    
    #set stratospheric water vapor 
    #(above 150 hPa set to the toga profile, so it isn't too moist)
    q_humid = q_humid*1e3
    q_humid[-3:] = q_toga[-3:]
        
    
    snd_out[:,0] = zlevs
    snd_out[:,1] = plevs/100.
    snd_out[:,2] = thetas
    snd_out[:,3] = q
    snd_out[:,4] = u
    snd_out[:,5] = v
    
    
    day0=0
    
    fname_out = 'snd_{0}RCE'.format(T_s)
    fname_out = 'snd_DRYRCE'
    head='z[m] p[mb] T[K] q[g/kg] u[m/s] v[m/s]'
    
    plt.figure(1)

    plt.plot(RH_new*100, plevs/1e2, 'r', label = 'dry')
    #plt.plot(RH_toga*100, plevs/1e2, 'b', label = 'humid')
    plt.xlabel('RH (%)')
    plt.ylabel('p (hPa)')
    ax = plt.gcf().gca()
    ax.set_yscale('log')
    plt.yticks([1000, 500, 250, 100, 50, 20])
    ax.set_ylim(50, 1000)
    ax.invert_yaxis()
    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.show()
    
    plt.figure(3)
    plt.plot(T, plevs/1e2)
    plt.xlabel('T (K)')
    plt.ylabel('p (hPa)')
    ax = plt.gcf().gca()
    ax.set_yscale('log')
    plt.yticks([1000, 500, 250, 100, 50, 20])
    ax.set_ylim(50, 1000)
    ax.invert_yaxis()
    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.show()
    
    plt.figure(5)
    plt.plot(q, plevs/1e2, 'r', label='dry')
    plt.plot(q_humid, plevs/1e2, 'b', label='humid')
    plt.xlabel(r'$q_v$ (g/kg)')
    plt.ylabel('p (hPa)')
    ax = plt.gcf().gca()
    ax.set_yscale('log')
    plt.yticks([1000, 500, 250, 100, 50, 20])
    ax.set_ylim(50, 1000)
    ax.invert_yaxis()
    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.legend()
    plt.show()
  
    
    np.savetxt(fpath + fname_out, snd_out, delimiter=' ', header=head, fmt='%3.4f')
    



    
    