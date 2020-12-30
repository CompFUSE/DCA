# Author: Giovanni Balduzzi
# Usage: python3 visualize_configuration.py <configuration file> <time stamp>

import numpy as np
import h5py
from sys import argv

if(len(argv) != 3) :
    print("Arguments: <file name> <time stamp>")
    exit(-1)

filename = argv[1]
stamp = argv[2]

file = h5py.File(filename,'r')[stamp]

times = file['times'][:]
sites = file['sites'][:, 0]
spins = file['hs_spin'][:]

# remove non-interacting spins
times = times[spins != 0]
sites = sites[spins != 0]
spins = spins[spins != 0]
