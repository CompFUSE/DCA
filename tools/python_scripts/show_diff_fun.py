# Prints the name of all groups and datasets stored in an HDF5 file.
#
# Usage: python show_hdf5.py <HDF5-output-file>
#
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)

import h5py
import sys
import numpy as np
import math
# Helper function
def print_name(name):
    print(name)

def compute_error(arr1, arr2):
	l1 = 0.
	l2 = 0.
	linf = 0.
	l1_error = 0.
	l2_error = 0.
	linf_error = 0.


	for i in range(len(arr2)):
		f_abs = abs(arr1[i])
		l1 += f_abs
		l2 += f_abs * f_abs
		if(linf<f_abs):
			linf=f_abs
			index_max_linf = i
		err = abs(arr1[i] - arr2[i])
		l1_error += err
		l2_error += err * err
		if(linf_error < err):
			linf_error = err
			index_max_linf_error = i
		if(arr1[i] == arr2[i]):
			if(arr1[i] != 0):
				print("non-zero equal at index: ", i, " val: ", arr1[i])
			else:
				print("zero equal at index: ", i, " val: ", arr1[i])
		elif (arr1[i] - arr2[i] < 5E-07):
				print("acceptable at index: ", i, "val: ", arr1[i])
		else:
			print("diff at index: ", i, "G4_1: ", arr1[i], " G4_2: ", arr2[i].real)
	
	l1_error /= l1
	l2_error = math.sqrt(l2_error / l2)
	linf_error /= linf	
	
	print("l1_error: ", l1_error)
	print("l2_error: ", l2_error)
	print("linf_error: ", linf_error)
	#print("index_max_linf: ", index_max_linf, "G4_1_f: ", G4_1_f[index_max_linf], " G4_2_f: ", G4_2_f[index_max_linf])
	#print("index_max_linf_error: ", index_max_linf_error, "G4_1_f: ", G4_1_f[index_max_linf_error], " G4_2_f: ", G4_2_f[index_max_linf_error])
	
# Get the filename from the command line.
filename1 = sys.argv[1]
filename2 = sys.argv[2]

# Open the file.
data1 = h5py.File(filename1,'r')
data2 = h5py.File(filename2,'r')

print("comparing file ", filename1, " and ", filename2)

# get G4
G4_1 = data1["functions/G4_PARTICLE_PARTICLE_UP_DOWN"][:]
G4_2 = data2["functions/G4_PARTICLE_PARTICLE_UP_DOWN"][:]

print("G4_1 shape: ", G4_1.shape)
print("G4_1 size: ", G4_1.size)
print("G4_2 shape: ", G4_2.shape)
print("G4_2 size: ", G4_2.size)

# compare apple to apple
assert(G4_1.shape == G4_2.shape)
assert(G4_1.size == G4_2.size)

np.set_printoptions(suppress=True)
# flatten n-d array to 1d
G4_1_f = G4_1.flatten();
G4_2_f = G4_2.flatten();


print("error in real part: ")
compute_error(G4_1_f.real, G4_2_f.real)
print("\n\nerror in imag part: \n\n")
compute_error(G4_1_f.imag, G4_2_f.imag)
