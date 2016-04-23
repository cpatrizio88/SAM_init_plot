import os
import glob

# converts all stat files in the current directory to netcdf files 
fpath = '/Users/cpatrizio/SAM6.10.8/OUT_STAT/'
binFile_names = glob.glob(fpath + '*.stat')
#path to bin3D2nc script, which converts a single bin3D file to .nc file
fn = fpath + 'stat2nc'

for f in binFile_names:
    cmd = fn + ' ' + f
    os.system(cmd)