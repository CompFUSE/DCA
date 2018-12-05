#!/bin/bash
#
# Loads all modules that are required to build DCA++ on ORNL's Summit supercomputer.
# A reset is done at the beginning to restore to the default programming environment on Summit.
#
# Usage: source summit_load_modules.sh

module reset
module load gcc/6.4.0
module load cuda
module load hdf5
module load fftw
module load cmake 
module load magma
module load netlib-lapack
module load essl

export CC=mpicc
export CXX=mpicxx
#export BLAS=$OLCF_ESSL_ROOT
#export LD_LIBRARY_PATH=$BLAS/lib:$LD_LIBRARY_PATH

