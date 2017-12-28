# Prints the name of all groups and datasets stored in an HDF5 file.
#
# Usage: python show_hdf5.py <HDF5-output-file>
#
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)

import h5py
import sys

# Helper function
def print_name(name):
    print(name)

# Get the filename from the command line.
filename = sys.argv[1]

# Open the file.
data = h5py.File(filename,'r')

# Recursively visit all objects (groups or datasets) and print their name.
data.visit(print_name)
