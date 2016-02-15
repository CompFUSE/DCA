# Usage: python plot.diff_self-energy.py input1.hdf5 input2.hdf5 ...

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
s1 = 0
s2 = 0
wn_cutoff = 32


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
wn = data_0['domains/frequency-domain/elements']


# Plot
for K_ind, K_vec in enumerate(K_cluster):
    plt.figure(K_ind)

    # Real part
    plt.subplot(211)
    for f_ind, filename in enumerate(filenames):

        data = h5py.File(filename,'r')
        self_energy = data['functions/Self_Energy/data']

        plt.plot(wn[len(wn)/2-wn_cutoff:len(wn)/2+wn_cutoff], self_energy[len(wn)/2-wn_cutoff:len(wn)/2+wn_cutoff, K_ind, s2, b2, s1, b1, 0],
                 marker=markers[f_ind], color=colors[f_ind], label=filename)

    # plt.xlabel('$\omega_n$')
    plt.ylabel(r'$\mathrm{Re} \, \Sigma(\vec{K}, i \omega_n)$')
    plt.legend(loc='upper right')
    plt.title(r'$\vec{K}=('+str(K_vec[0])+', '+str(K_vec[1])+')$')
        
    # Imaginary part
    plt.subplot(212)
    for f_ind, filename in enumerate(filenames):
        
        data = h5py.File(filename,'r')
        self_energy = data['functions/Self_Energy/data']

        plt.plot(wn[len(wn)/2-wn_cutoff:len(wn)/2+wn_cutoff], self_energy[len(wn)/2-wn_cutoff:len(wn)/2+wn_cutoff, K_ind, s2, b2, s1, b1, 1],
                 marker=markers[f_ind], color=colors[f_ind], label=filename)
        
    plt.xlabel('$\omega_n$')
    plt.ylabel(r'$\mathrm{Im} \, \Sigma(\vec{K}, i \omega_n)$')
    plt.legend(loc='upper right')
    plt.savefig('self-energy_K='+str(K_ind)+'.pdf')
