# Plots the eigenvector of the leading eigenvalue vs. fermionic Matsubara frequency for each
# cluster (DCA) or tp-host (DCA+) momentum.
#
# Usage: python plot_eigenvector.py <analysis-HDF5-output-file>
#
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)

from matplotlib import pyplot as plt
import h5py
import os
import sys

# Get the filename from the command line.
filename = sys.argv[1]

# Open the file.
data = h5py.File(filename,'r')

# Fermionic Matsubara frequency mesh
wn = data['domains/vertex-frequency-domain (COMPACT)/elements'] [:]
# print(wn)

# Eigenvector momentum grid
k = data['domains/LATTICE_TP/MOMENTUM_SPACE/elements/data'][:]
# print(k)

# Eigenvectors of leading eigenvalues
# shape: [wn, k, b2, b1, ev_index, real/imag],
#        where ev_index labels the eigenvectors of the 10 leading eigenvalues.
#        ev_index = 0 corresponds to the largest eigenvalue.
eigenvecs = data['analysis-functions/leading-eigenvectors/data'][:]
# print(eigenvecs.shape)

# Parameters for plot
# Band and eigenvalue indices
b1 = 0
b2 = 0
ev_index = 0

plt.figure()

for k_ind, k_vec in enumerate(k):
    # Plot real part.
    plt.plot(wn, eigenvecs[:, k_ind, b2, b1, ev_index, 0], 'o-', ms=2,
             label=r'$\mathbf{K}$'+'=('+'{0:.3g}'.format(k_vec[0])+', '+'{0:.3g}'.format(k_vec[1])+')')

plt.legend(loc='best')
plt.xlabel('$\omega_n$')
plt.ylabel(r'$\phi(\mathbf{K}, \omega_n)$')
plt.savefig('leading_eigenvector.pdf')
