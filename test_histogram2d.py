import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm

r = np.linspace(0, 500e3, 100)
z = np.linspace(0, 25e3, 64)

rr, zz = np.meshgrid(r, z)

weights = rr + zz*50

nbins=(100, 10)

fieldsums, redges, zedges = np.histogram2d(rr.ravel(), zz.ravel(), bins=nbins, weights=weights.ravel())

nr, redges, zedges = np.histogram2d(rr.ravel(), zz.ravel(), bins=nbins)

rbin_centers = (redges[1:] + redges[:-1])/2.
zbin_centers = (zedges[1:] + zedges[:-1])/2.
    
rr, zz = np.meshgrid(rbin_centers, zbin_centers)

fieldmeans = fieldsums/nr

plt.figure()
plt.pcolormesh(rr/(1e3), zz/(1e3), np.transpose(fieldmeans), cmap=cm.RdYlBu_r)
plt.show()