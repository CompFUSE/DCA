# Usage: python plot.diff_sp-greens-function.py input1.hdf5 input2.hdf5 ...

import h5py
import sys

import numpy as np
import scipy as sp

from matplotlib import pyplot as plt
from matplotlib import colors


def printname(name):
    print name

    
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
markers = ['o', 's', '^', 'v', '<', '>', '*']


# Parameters
b1 = 0
b2 = 0
b3 = 0
b4 = 0


filenames = sys.argv[1:]
# print filenames

if len(filenames) == 0:
    print "No input files."
    sys.exit(0)

if len(filenames) > 7:
    print "Too many input files."
    sys.exit(0)

    
# Read cluster K-points and Matsubara frequencies from first input file
data_0 = h5py.File(filenames[0],'r')
# data_0.visit(printname)

K_cluster = data_0['domains/CLUSTER/MOMENTUM_SPACE/elements/data']
wn = data_0['domains/vertex-frequency-domain (COMPACT)/elements']


# Plot
for K_ind, K_vec in enumerate(K_cluster):
    plt.figure(K_ind)

    # Real part
    plt.subplot(211)
    for f_ind, filename in enumerate(filenames):

        data = h5py.File(filename,'r')
        G4 = data['functions/G4_k_k_w_w/data']

        G4_Re = []
        for w_ind, w in enumerate(wn):
            G4_Re.append(G4[w_ind, w_ind, K_ind, K_ind, b4, b3, b2, b1, 0])
            
        plt.plot(wn, G4_Re, marker=markers[f_ind], color=colors[f_ind], label=filename)

    # plt.xlabel('$\omega_n$')
    plt.ylabel(r'$\mathrm{Re} \, G4(\vec{K}, \vec{K}, i \omega_n, i \omega_n)$')
    plt.legend(loc='upper right')
    plt.title(r'$\vec{K}=('+str(K_vec[0])+', '+str(K_vec[1])+')$')
        
    # Imaginary part
    plt.subplot(212)
    for f_ind, filename in enumerate(filenames):
        
        data = h5py.File(filename,'r')
        G4 = data['functions/G4_k_k_w_w/data']

        G4_Im = []
        for w_ind, w in enumerate(wn):
            G4_Im.append(G4[w_ind, w_ind, K_ind, K_ind, b4, b3, b2, b1, 1])
            
        plt.plot(wn, G4_Im, marker=markers[f_ind], color=colors[f_ind], label=filename)
        
    plt.xlabel('$\omega_n$')
    plt.ylabel(r'$\mathrm{Im} \, G4(\vec{K}, \vec{K}, i \omega_n, i \omega_n)$')
    plt.legend(loc='upper right')
    plt.savefig('tp-greens-function_K='+str(K_ind)+'.pdf')
