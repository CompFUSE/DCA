#!/bin/bash
#
# Loads all modules that are required to build DCA++ on ORNL's Summit supercomputer.
# A reset is done at the beginning to restore to the default programming environment on Summit.
#
# Usage: source summit_load_modules.sh

module reset
module load gcc/8.1.1
module load cuda
module load magma/2.5.1
module load hdf5
module load fftw
module load cmake
module load netlib-lapack
module load essl

export CC=mpicc
export CXX=mpicxx
