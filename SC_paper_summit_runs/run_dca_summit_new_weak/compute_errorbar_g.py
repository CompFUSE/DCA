#/usr/bin/python3

# Usage: python compute_errorbar.py <HDF5-output-file>

# Author: Giovanni Balduzzi

import h5py
import sys
import numpy as np

if len(sys.argv) < 2:
    print("Usage: python compute_errorbar.py <HDF5-output-file>")
    exit

filename = sys.argv[1]
print("Reading from ", filename)

file = h5py.File(filename,'r')
print("Loading error.")
G_err = file["functions/G_k_w-error/data"][:]
print("Loading value.")
G = file["functions/cluster_greens_function_G_k_w/data"][:]
assert(len(G) == len(G_err))

print("Computing norms")

def norms(f):
    l1 = np.sum(f)
    l2 = np.sum(f * f)
    lmax = np.amax(f)

    return np.array([l1, np.sqrt(l2), lmax])

err = norms(G_err)
norm = norms(G)
#ratio = norms(G4_err/G4)
 

print("Error norm: l1 l2 l_inf")
print(err)
print("G norm: l1 l2 l_inf")
print(norm)
print("ratio of norms: l1 l2 l_inf")
print(err/norm)
#print("norm of ratios: l1 l2 l_inf")
#print(ratio)

