import h5py
import sys
import numpy as np


ns = [100, 500, 1000, 2000, 4608]

def norm(f):
    l2  = np.sum(f * f)
    return l2

data = []
for n in ns :
    file = h5py.File('nodes_' + str(n) + '/dca.hdf5', 'r')
    g_err = norm(file["functions/G_k_w-error/data"][:])
    g4_err =norm(file["functions/G4-error/data"][:])

    data.append([n, g_err, g4_err])

data = np.array(data)
np.savetxt('errors.txt', data, header='n\tg_err_l2\tg4_err_l2')
