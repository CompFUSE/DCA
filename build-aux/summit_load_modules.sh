#!/bin/bash
#
# Loads all modules that are required to build DCA++ on ORNL's Summit supercomputer.
# A reset is done at the beginning to restore to the default programming environment on Summit.
#
# Usage: source summit_load_modules.sh

module reset
module load gcc/12.1.0
module load cuda/11.5.2 # ldd shows magma is built with this cuda
module load magma/2.6.2
module load hdf5
module load fftw
module load cmake/3.21.3
module load netlib-lapack
module load essl
module load adios2
module load ninja

export CC=mpicc
export CXX=mpic++
