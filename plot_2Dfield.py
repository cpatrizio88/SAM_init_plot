import numpy as np
from netCDF4 import Dataset
import site
site.addsitedir('/Users/cpatrizio/Dropbox/research/code/thermlib/')
site.addsitedir('/Users/cpatrizio/Dropbox/research/')
from findT import *
from theta import *
import matplotlib.pyplot as plt
import glob

fpath =  '/Users/cpatrizio/SAM6.10.8/OUT_2D_BACKUP/'
nc_in = glob.glob(fpath + '*.nc')[0]

ncfile_in = Dataset(nc_in)

var = ncfile_in.variables

x = var['x'][:]
y = var['y'][:]
t = var['time'][:]



    
    
    



