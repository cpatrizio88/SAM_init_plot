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
    
def radprof(field, xx, yy, center, nbins):
    """
    compute radial profile from the 'center' point for 2D field
    assumes square, periodic domain, given by x,y.
    
    """
    
    width = xx[0,-1]
    delx = np.abs(xx - center[0])
    dely = np.abs(yy - center[1])
    delx[delx > width/2] = width - delx[delx > width/2]
    dely[dely > width/2] = width - dely[dely > width/2]
    r = np.sqrt(np.power(delx, 2) + np.power(dely, 2))
    #r = r.astype(np.int)
    
    #min_r = np.min(r)
    
    
    fieldsums, redges = np.histogram(r, bins=nbins, weights=field)
    nr, redges = np.histogram(r, bins=nbins)
    
   # nr = np.bincount(r.ravel())
   # radialsums = np.bincount(r.ravel(), weights=field.ravel())
    
    fieldmeans = fieldsums/nr

    #r_max = np.max(r)
    #radialsum, radius = np.histogram(r, weights=field, bins=r_max/3)
    
    return (redges, fieldmeans)
    
def radprof3D(field, xx, yy, z, center, nbins):
    """
    compute radial profile from the 'center' point for 3D field
    assumes square, periodic domain, given by x,y.
    
    """
    #nx = len(xx[:,0])
    #ny = len(yy[:,0])
    nz = z.size
    width = xx[0,-1]
    delx = np.abs(xx - center[0])
    dely = np.abs(yy - center[1])
    delx[delx > width/2] = width - delx[delx > width/2]
    dely[dely > width/2] = width - dely[dely > width/2]
    r = np.sqrt(np.power(delx, 2) + np.power(dely, 2))
    
    r = r.ravel()
    field = field.reshape(nz, r.size)
    #r = np.repeat(r[np.newaxis,:], z.size, axis=0)
    rr, zz = np.meshgrid(r, z)
  
    fieldsums, redges, zedges = np.histogram2d(rr.ravel(), zz.ravel(), bins=nbins, weights=field.ravel())
    nr, redges, zedges = np.histogram2d(rr.ravel(), zz.ravel(), bins=nbins)
    
    fieldmeans = fieldsums/nr
    
    return (redges, zedges, fieldmeans)
    
def fracclusterarea(name, varis, nave, t, a=1, W500crit=0.025):
    """
    computes the fractional area of a convective region
    if name == 'PW', 
       convective regions are defined as points where PW > mean(PW) + a*std(PW).
       assume dx=dy
    if name == 'W500' 
       convective regions are defined as points where the time-mean vertical velocity is greater
       than zero
       
    nave is the number of days to average over
    
    a is a threshold factor 
    

    """
    tvari = varis['PW']
    if len(tvari.shape) > 2:
        nx = tvari.shape[1]
        ny = tvari.shape[2]
        totpoints = nx*ny
    else:
        nx = tvari.shape[1]
        totpoints = nx
    #
    #if (name == 'MIX'):
    #    P = varis['Prec'][:]
    #    W500 = varis['W500'][:]
    #    Pavefield = np.mean(P[t-nave:t,:,:], axis=0)
    #    W500avefield = np.mean(W500[t-nave:t,:,:], axis=0)
    #    Pbar = np.mean(Pavefield)
    #    Pstd = np.mean(Pavefield)
    #    Pcrit = Pbar + a*Pstd
    #    convecpoints = np.bitwise_and(Pavefield >= Pcrit, W500avefield >= W500crit)
    #    count = len(Pavefield[convecpoints])
    #    return count/float(totpoints)
    
    if (name == 'PrecNEW'):
        P = varis['Prec'][:]
        Pavefield = np.mean(P[t-nave:t,:,:], axis=0)
        Pbar = np.mean(Pavefield)
        Phat = np.mean(Pavefield[Pavefield > 0])
        return Pbar/Phat
        
    vari = varis[name][:]
    avefield = np.mean(vari[t-nave:t,:,:], axis=0)

    
    if (name == 'PW'):
        PWbar = np.mean(avefield)
        PWstd = np.std(avefield)
        count = len(avefield[avefield >= PWbar + a*PWstd])
        return count/float(totpoints)
    if (name == 'Prec'):
        Pbar = np.mean(avefield)
        Pstd = np.std(avefield)
        count = len(avefield[avefield >= Pbar + a*Pstd])
        return count/float(totpoints)
    #NOT WORKING ?
    #if (name == 'W500'):
    #    print 'W500crit', W500crit
    #    print 'totpoints', totpoints
    #    count = len(avefield[avefield >= W500crit])
    #    return count/float(totpoints)
        
    
    return ()
        
        
    
    
    
    
    
    

    