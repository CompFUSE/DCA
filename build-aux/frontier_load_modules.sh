#!/bin/bash
#
# Loads all modules that are required to build DCA++ on ORNL's Frontier.
# A reset is done at the beginning to restore to the default programming environment on Frontier.
# This is for development only at this point.
#
# Usage: source frontier_load_modules.sh


module reset
module load gcc/11.2.0
module load openblas
module load hdf5
module load amd-mixed/5.3.0
module load fftw
module load cmake
module load ninja
module use /sw/frontier/magma/lmod/Core/cce-15.0.0
module load magma/2.7.2
module load craype-accel-amd-gfx90a

export CC=mpicc
export CXX=mpicxx

export OPENBLAS_ROOT=${OLCF_OPENBLAS_ROOT}
export HDF5_ROOT=${OLCF_HDF5_ROOT}
export MAGMA_ROOT=/sw/frontier/magma/opt/linux-sles15-zen3/cce-15.0.0/magma-2.7.2-x7o7sph6npwb73t2leetbgpbypqyhtz6
export FFTW_PATH=${OLCF_FFTW_ROOT}
export LD_PRELOAD=/opt/cray/pe/lib64/cce/libtcmalloc_minimal.so.1
