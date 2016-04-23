import numpy as np
import site
site.addsitedir('/Users/cpatrizio/Dropbox/research/code/thermlib/')
site.addsitedir('/Users/cpatrizio/Dropbox/research/')
from findT import *


fpath_in =  '/Users/cpatrizio/SAM6.10.8/TOGA/'
fname_in = 'grd'


grd_in = np.loadtxt(fpath_in + fname_in, skiprows=0)

#64 vertical levels
zlevs = grd_in[:,0]

#32 vertical levels
zlevs_32 = zlevs[1::2]

fpath_out = '/Users/cpatrizio/SAM6.10.8/RCE_IDEAL/'
fname_out = 'grd_32'

np.savetxt(fpath_out + fname_out, zlevs_32, fmt='%3.4f')