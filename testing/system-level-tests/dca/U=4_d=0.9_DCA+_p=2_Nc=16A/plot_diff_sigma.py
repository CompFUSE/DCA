
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

filename = "./T=0.25/data.DCA+_sp.hdf5"

# Parameters
b1 = 0
b2 = 0
s1 = 0
s2 = 0
wn_cutoff = 32

mode = "DCA+"
Nc = "16A"
period = "2"
    
# Read cluster K-points and Matsubara frequencies from first input file
# find K=(pi,0)

data = h5py.File(filename,'r')

wn    = data["domains"]["frequency-domain"]["elements"][240:272]  # 240:272 chosen from 512 sp frequencies in sp hdf5 file
k_dmn = data["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"]["data"]

shape = k_dmn.shape

Kxsp = k_dmn[:,0]  
Kysp = k_dmn[:,1]  

for iK in range(0,shape[0]):
    if abs(Kysp[iK]<1.e-3) and abs(Kxsp[iK]-pi)<1.e-3:
       iKpi0 = iK

# Plot reference data first
a = loadtxt("./sigma_w_K=pi0_reference")
plt.plot(a[:,0], a[:,1], 'kp-', label="sigma_w_K=pi0_reference")
        
# Plot generated Im sigma(K=(pi,0),iwn)
sigma = data['functions/Self_Energy/data']
plt.plot(wn, sigma[240:272, iKpi0, s2, b2, s1, b1, 1], 'bo-', label=filename)

plt.xlabel('$\omega_n$')
plt.ylabel(r'$\mathrm{Im} \, \Sigma(\vec{K}, i \omega_n)$')
plt.legend(loc='upper right')
plt.savefig('check_self-energy_K=pi0_'+mode+'_p='+period+'_Nc='+Nc+'.pdf')

