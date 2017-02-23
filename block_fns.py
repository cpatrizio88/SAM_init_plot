import numpy as np
import collections
from thermolib.constants import constants

c = constants()

def blockave1D(field, db):
    """
    Takes a 1D field with nx points and divides into blocks. 
    The new field has length nx/db.   
    
    """
    
    nx = field.shape[0]
    nxblock = nx // db
    
    blockfield = np.average(np.split(field, nxblock, axis=0), axis=1)
    
    return blockfield


def blockave2D(field, db):
    """
    Takes a nx x ny field and divides the field into blocks - the new field has
    dimensions nx' x ny' where nx' = nx/db, ny' = ny/db
    db is the block width in number of grid cells. 
    the field is averaged over the block area (db points) 
    
    in the case of 3D field, averaging is only performed in horizontal direction. 
    
    assumes nx = ny = even integer.
    db must be a multiple of nx
    
    
    """
    
    nx = field.shape[0]
    ny = field.shape[1]
    
    nxblock = nx // db
    nyblock = ny // db
    
    #split up field column-wise, take average row-wise. then split up resulting field row-wise, and take average column-wise.
    blockfield = np.average(np.split(np.average(np.split(field, nxblock, axis=1), axis=-1), nyblock, axis=1), axis=-1)
    
    return blockfield
    
    
def blockave3D(field, db):
  
    
    nz = field.shape[0]
    nx = field.shape[1]
    ny = field.shape[2]
    
    nxblock = nx // db
    nyblock = ny // db
    
    blockfield = np.zeros((nz, nxblock, nyblock))
    
    for i in range(nz):
        field_z = field[i,:,:]
        blockfield[i,:,:] = np.average(np.split(np.average(np.split(field_z, nxblock, axis=1), axis=-1), nyblock, axis=1), axis=-1)
    
    return blockfield
    
def blocksum3D(field, db):
  
    
    nz = field.shape[0]
    nx = field.shape[1]
    ny = field.shape[2]
    
    nxblock = nx // db
    nyblock = ny // db
    
    blockfield = np.zeros((nz, nxblock, nyblock))
    
    for i in range(nz):
        field_z = field[i,:,:]
        blockfield[i,:,:] = np.sum(np.split(np.sum(np.split(field_z, nxblock, axis=1), axis=-1), nyblock, axis=1), axis=-1)
    
    return blockfield
    
    
    
def blocksort1D(sfield, ofield, db):
    
    """
    Takes two 1D fields of length nx, and sorts according to sfield.
    
    """
    
    nx = sfield.shape[0]
    nxblock = nx // db
    
    blocksfield = np.average(np.split(sfield, nxblock, axis=0), axis=1)
    blockofield =  np.average(np.split(ofield, nxblock, axis=0), axis=1)
    
    d = dict(zip(blocksfield, blockofield))
    od = collections.OrderedDict(sorted(d.items()))
    
    return od
    
def blocksort2D_1Din(sfield, ofield, db):
    
    """
    Takes a 1D field of length nx, sfield, and 2D (x-z) of length nz by nx and sorts according to sfield (preserves z shape).
    
    """
    
    nx = sfield.shape[0]
    nz = ofield.shape[0]
    nxblock = nx // db
    
    blocksfield = np.average(np.split(sfield, nxblock, axis=0), axis=1)
    
    blockofield = np.zeros((nz, nxblock))
    
    for i in range(nz):
        field_z = ofield[i,:]
        blockofield[i,:] = np.average(np.split(field_z, nxblock, axis=0), axis=1)
        
    sorter = blocksfield.argsort()
    blocksfield = np.sort(blocksfield)
    blockofield = blockofield[:,sorter]
    
    return (blocksfield, blockofield)
    
    
def blocksort2D(sfield, ofield, db):
    """
    Takes two nx x ny fields and divides them into blocks - the new fields have
    dimensions nx' x ny' where nx' = nx/db, ny' = ny/db
    db is half the block width in number of grid cells. 
    the fields are averaged over the block area (db points) and then 
    ofield is sorted according to sfield (spatial structure is lost)
    
    the returned value is a dictionary with sfield as the key and ofield as the value 
    
    assumes nx = ny = even integer.
    db must be a multiple of nx
    
    """
    nx = sfield.shape[0]
    ny = sfield.shape[1]
    
    nxblock = nx / db
    nyblock = ny / db
    
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
    
def blocksort3D(sfield, ofield, db):
    """
    Takes a nx x ny field and a nz x nx x ny field and divides them into horizontal blocks - the new fields have
    dimensions nx' x ny' where nx' = nx/db, ny' = ny/db
    db is half the block width in number of grid cells. 
    the fields are averaged over the horizontal block area (db points) and then 
    ofield is sorted according to sfield 
    (horizontal spatial structure is lost, only retain vertical structure.. so the resulting ofield is sorted profiles)
    
    the returned value is a tuple with the sorted fields (sfield, ofield)
     
    assumes nx = ny = even integer.
    db must be a multiple of nx
    
    """
    
    nz = ofield.shape[0]
    nx = ofield.shape[1]
    ny = ofield.shape[2]
    
    nxblock = nx // db
    nyblock = ny // db
    
    blocksfield = np.average(np.split(np.average(np.split(sfield, nxblock, axis=1), axis=-1), nyblock, axis=1), axis=-1)
    
    blockofield = np.zeros((nz, nxblock, nyblock))
    
    for i in range(nz):
        field_z = ofield[i,:,:]
        blockofield[i,:,:] = np.average(np.split(np.average(np.split(field_z, nxblock, axis=1), axis=-1), nyblock, axis=1), axis=-1)
        
    blocksfield = blocksfield.flatten()
    
    blockofield = blockofield.reshape((nz, nxblock*nyblock))
    
    sorter = blocksfield.argsort()
    blocksfield = np.sort(blocksfield)
    blockofield = blockofield[:,sorter]
    
    return (blocksfield, blockofield)
    
def blockxysort2D(sfield, xx, yy, db):
    """
    Takes a nx x ny field and divides it into blocks - the new field has
    dimensions nx' x ny' where nx' = nx/db, ny' = ny/db
    db is half the block width in number of grid cells. 
    the field is averaged over the block area (db points) and then 
    ofield is sorted according to sfield along with the (x,y) coords of the block
    
    the returned value is a dictionary with sfield as the key and (x,y) as value
    
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
    xblock = np.average(np.split(np.average(np.split(xx, nxblock, axis=1), axis=-1), nyblock, axis=1), axis=-1)
    yblock = np.average(np.split(np.average(np.split(yy, nxblock, axis=1), axis=-1), nyblock, axis=1), axis=-1)
    
    xblock = xblock.flatten()
    yblock = yblock.flatten()
    xyblock = zip(xblock, yblock)
    blocksfield = blocksfield.flatten()

    d = dict(zip(blocksfield, xyblock))
    
    od = collections.OrderedDict(sorted(d.items()))
    
    return od
    
def xysort(sfield, x, y):
    x = x.flatten()
    y = y.flatten()
    sfield = sfield.flatten()
    xy = zip(x, y)
    d = dict(zip(sfield, xy))
    od = collections.OrderedDict(sorted(d.items()))
    return od
    
def vertint(field, p):
    """
    vertically integrates a 3D field over pressure. 
    if hydrostatic atmoshpere, this is equivalent to density weighted vertical integral.
    ignore the last level of field when integrating to have consistent shape with diff(p)
    """
    nz = field.shape[0]
    nx = field.shape[1]
    ny = field.shape[2]
    dp3D= np.ones((nz-1, nx, ny))
    dp3D= np.transpose(dp3D)
    dp = -np.diff(p)
    dp3D[:,:,:]=dp
    dp3D = np.transpose(dp3D)
    
    fieldhat = (1/c.g)*np.sum(np.multiply(field[:-1,:,:], dp3D), axis=0)
    return fieldhat
    
def vertint2D(field, p):
    """
    vertically integrates a 2D (x-z) field over pressure. 
    if hydrostatic atmoshpere, this is equivalent to density weighted vertical integral.
    ignore the last level of field when integrating to have consistent shape with diff(p)
    """
    nz = field.shape[0]
    nx = field.shape[1]
    dp2D= np.ones((nz-1, nx))
    dp2D= np.transpose(dp2D)
    dp = -np.diff(p)
    dp2D[:,:]=dp
    dp2D = np.transpose(dp2D)
    
    fieldhat = (1/c.g)*np.sum(np.multiply(field[:-1,:], dp2D), axis=0)
    return fieldhat
    
    

  
    
    
    
    
    
        
        
    
    
    
    