# Computes the transition temperature Tc from the temperature dependence of the leading
# Bethe-Salpeter eigenvalue.
#
# Usage: python compute_tc.py T=*
#
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)

import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt
import h5py
import os
import sys


################################################################################

# Computes the temperature at which an instability occurs, i.e. the temperature T where the leading
# eigenvalue eigval crosses 1.
# Uses a fit function of the form eigval(T) = p0/(T-p1)^p2.
# The transition temperature Tc is then given by Tc = p1^(1/p0) + p1.
def computeTransitionTemp(T, eigval):
    print('\nTemperature/eigenvalue pairs for fit:')
    for T_ind, T_val in enumerate(T):
        print(str(T_val) + '\t' + str(eigval[T_ind]))

    fitfunc = lambda p, x: p[0] / pow((x-p[1]), p[2])  # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y        # Distance to the target function

    p0 = [1., 0., 1.]  # Initial guess for the parameters
    p, success = optimize.leastsq(errfunc, p0[:], args=(T, eigval))

    Tc = pow(p[0], 1./p[2]) + p[1]
    print('\nTc = ' + '{0:.3g}'.format(Tc))

    T_fine = np.linspace(T[0], T[-1], 100)
    l_fine = fitfunc(p, T_fine)

    return Tc, T_fine, l_fine


################################################################################

dirs = sys.argv[1:]  # T=... directories

T = []
eigval = []

# Read leading eigenvalue for each temperature.
for d in dirs:
    filename = d + '/analysis.hdf5'

    if (os.path.isfile(filename)):
        T.append(float(d[2:]))

        print('Reading ' + filename)
        data = h5py.File(filename,'r')

        # Store real part of leading eigenvalue (imaginary part = 0).
        # Eigenvalues are sorted w.r.t. size in decreasing order.
        leading_eigenvalues = data['analysis-functions']['leading-eigenvalues']['data'][:]
        eigval.append(leading_eigenvalues[0][0])

        data.close()

# Compute the transition temperature Tc.
Tc, T_fine, eigval_fine = computeTransitionTemp(T, eigval)

# Plot temperature dependence of leading eigenvalue.
filename = 'eigval_vs_temp.pdf'
print('\nPlotting temperature dependence of leading eigenvalue: ' + filename)

xmin = T_fine[0]-0.005
xmax = T_fine[-1]+0.005

plt.plot(T_fine, eigval_fine, '--', label=r'$T_c$ = '+'{0:.3g}'.format(Tc))
plt.plot(T, eigval, 'o')
plt.hlines(1., xmin, xmax, 'k')
plt.xlim(xmin, xmax)
plt.xticks([0.07, 0.08, 0.09, 0.1])
plt.xlabel(r'$T/t$')
plt.ylabel(r'$\lambda_d$')
plt.legend(loc='best')
plt.grid()
plt.savefig(filename)
