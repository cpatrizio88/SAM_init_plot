import numpy as np
import site
#site.addsitedir('\\Users\\Casey\\Dropbox\\research\\code\\thermlib')
site.addsitedir('/Users/cpatrizio/Dropbox/research/code/thermlib/')
from wsat import wsat


rho_w = 1000 #density of water
g=9.81 #gravitational acceleration
p_s = 1015


#function that takes netCDF variables, block-sorts the field given 
#by vname by CRH (column-relative humidity) at a given time i (the index of the time). 
#The blocks consist of nx grid points (in both x and y direction), 
#and the field is spatially averaged over the block, and daily averaged.

def CRH_blocksort(varis, vname, nx, ti, delt):
    x = varis['x'][:]
    y = varis['y'][:]
    z = varis['z'][:]
    p = varis['p'][:]
    p = np.concatenate(([p_s], p))
    p = p*1e2
    qv = varis['QV'][:]
    T = varis['TABS'][:]
    T_t = T[ti,:,:,:]
    qv_t = qv[ti,:,:,:]
    delp3D = np.ones((x.size, y.size, z.size))
    delp = -np.diff(p)
    delp3D[:,:,:] = delp
    delp3D = np.transpose(delp3D)
    
    PW = 1/(g*rho_w)*np.sum(np.multiply(delp3D, qv_t), axis=0)
    
    SPW = np.zeros(PW.shape)
    
    for i, plev in enumerate(p):
        print i
        SPW[:,:] = 1/(g*rho_w)*np.sum(np.multiply(plev, wsat(T_t[i], plev)))
        
    CRH = PW/SPW
    return CRH
    

    #SPW_fn = lambda p: (1./(g*rho_w)*wsat(findT(T_s, p), p))
    #PW = integrate.quad(SPW_fn, p_t, p_BL)[0] + (1./(g*rho_w))*wsat(T_BL, (p_BL + p_s)/2.)*p_BL






