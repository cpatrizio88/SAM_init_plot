import numpy as np
import matplotlib.pyplot as plt

import matplotlib

matplotlib.rcParams.update({'font.size': 26})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 2.5})
matplotlib.rcParams.update({'legend.fontsize': 22})


fout = '/Users/cpatrizio/Google Drive/figures/'


r = np.linspace(-1500, 1500.0, 1000)

r1 = 0
r2 = 300

K = 1e-5

A = 0.01

ur1 = -A*(r-r1)*np.exp(-K*(r-r1)**2)

ur2 = -A*(r-r2)*np.exp(-K*(r-r2)**2)

umerge = ur1 + ur2




plt.axvline(r1, color='r', alpha=0.1)
plt.axvline(r2, color='b', alpha=0.1)
plt.plot(r, ur1, linewidth=4, color = 'r', alpha=0.6, label = '$u_1$')
plt.plot(r, ur2, linewidth=4, color = 'b', alpha=0.6, label = '$u_2$')
plt.plot(r, umerge, color='k', label = r'$u_1 + u_2$' )
plt.axhline(0, color='k', alpha=0.1)
plt.text(r2-10, -2.12, r'$r_2$')
plt.xlabel(r'$r$ (km)')
plt.ylabel('surface wind radial velocity (m/s)')
plt.title(r'$r_1$ = {:1.0f} km, $r_2$ = {:3.0f} km'.format(r1,r2))
plt.legend(loc='best')
plt.xlim((0, 1500))
plt.savefig(fout + 'idealmerging_r2{:3.0f}km.pdf'.format(r2))
plt.close()

#plt.show()



