import numpy as np

def raddist(xx, yy, center):
    """
    compute radial distance from the 'center' point
    assumes square, periodic domain, given by x,y.
    
    """
    width = xx[0,-1]
    delx = np.abs(xx - center[0])
    dely = np.abs(yy - center[1])
    delx[delx > width/2] = width - delx[delx > width/2]
    dely[dely > width/2] = width - dely[dely > width/2]
    r = np.sqrt(np.power(delx, 2) + np.power(dely, 2))
    r = r.astype(np.int)
    return r
    
def radprof(field, xx, yy, center):
    """
    compute radial profile from the 'center' point
    assumes square, periodic domain, given by x,y.
    
    """
    
    width = xx[0,-1]
    delx = np.abs(xx - center[0])
    dely = np.abs(yy - center[1])
    delx[delx > width/2] = width - delx[delx > width/2]
    dely[dely > width/2] = width - dely[dely > width/2]
    r = np.sqrt(np.power(delx, 2) + np.power(dely, 2))
    r = r.astype(np.int)
    
    min_r = np.min(r)
    
    
    nr = np.bincount(r.ravel())
    radialsums = np.bincount(r.ravel(), weights=field.ravel())
    

    #r_max = np.max(r)
    #radialsum, radius = np.histogram(r, weights=field, bins=r_max/3)
    
    return (nr, radialsums)

    