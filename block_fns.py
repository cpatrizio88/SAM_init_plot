import numpy as np
import collections


def blockave2D(field, db):
    """
    Takes a nx x ny field and divides the field into blocks - the new field has
    dimensions nx' x ny' where nx' = nx/db, ny' = ny/db
    db is half the block width in number of grid cells. 
    the field is averaged over the block area (db points) 
    
    in the case of 3D field, averaging is only performed in horizontal direction. 
    
    assumes nx = ny = even integer.
    db must be a multiple of nx
    
    
    """
    
    nx = field.shape[0]
    ny = field.shape[1]
    
    nxblock = nx // db
    nyblock = ny // db
    
    #tave_field = np.mean(field[ti-ntave:ti,:,:])
    #tave_field = np.squeeze(tave_field)
    
    #split up field column-wise, take average row-wise. then split up resulting field row-wise, and take average column-wise.
    blockfield = np.average(np.split(np.average(np.split(field, nxblock, axis=1), axis=-1), nyblock, axis=1), axis=-1)
    
    return blockfield
    
    
def blockave3D(field, db):
  
    
    nz = field.shape[0]
    nx = field.shape[1]
    ny = field.shape[2]
    
    nxblock = nx // db
    nyblock = ny // db
    
    #tave_field = np.mean(field[ti-ntave:ti,:,:,:])
    #tave_field = np.squeeze(tave_field)
    
    blockfield = np.zeros((nz, nxblock, nyblock))
    
    for i in range(nz):
        field_z = field[i,:,:]
        blockfield[i,:,:] = np.average(np.split(np.average(np.split(field_z, nxblock, axis=1), axis=-1), nyblock, axis=1), axis=-1)
    
    return blockfield


def blocksort2D(sfield, ofield, db):
    """
    Takes two nx x ny fields and divides them into blocks - the new fields have
    dimensions nx' x ny' where nx' = nx/db, ny' = ny/db
    db is half the block width in number of grid cells. 
    the fields are averaged over the block area (db points) and then 
    ofield is sorted according to sfield (spatial structure is lost)
    
    the returned value is a dictionary with sfield as the key and ofield as the value 
    
    in the case of 3D field, averaging is only performed in horizontal direction 
    (so the value of the dictionary will be a 1D array). 
    
    assumes nx = ny = even integer.
    db must be a multiple of nx
    
    """
    nx = sfield.shape[0]
    ny = sfield.shape[1]
    
    nxblock = nx // db
    nyblock = ny // db
    
    #tave_field = np.mean(field[ti-ntave:ti,:,:])
    #tave_field = np.squeeze(tave_field)
    
    #split up field column-wise, take average row-wise. then split up resulting field row-wise, and take average column-wise.
    blocksfield = np.average(np.split(np.average(np.split(sfield, nxblock, axis=1), axis=-1), nyblock, axis=1), axis=-1)
    
    blockofield = np.average(np.split(np.average(np.split(ofield, nxblock, axis=1), axis=-1), nyblock, axis=1), axis=-1)
    
    blocksfield = blocksfield.flatten()
    blockofield = blockofield.flatten()
    
    d = dict(zip(blocksfield, blockofield))
    
    od = collections.OrderedDict(sorted(d.items()))
    
    return od
    
def vert_int(field, p):
    """
    vertically integrates a field over pressure. 
    if hydrostatic atmoshpere, this is equivalent to density weighted vertical integral.
    
    """
    dp3D= np.ones(field.shape)
    dp3D= np.transpose(dp3D)
    dp = -np.diff(p)
    dp3D[:,:,:]=dp
    dp3D = np.transpose(dp3D)
    
    fieldhat = np.sum(np.multiply(field, dp3D), axis=0)
    return fieldhat
    
    
    

  
    
    
    
    
    
        
        
    
    
    
    