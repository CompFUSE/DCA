#!/bin/bash
#
# Loads all modules that are required to build DCA++ on ORNL's Frontier.
# A reset is done at the beginning to restore to the default programming environment on Frontier.
# This is for development only at this point.
#
# Usage: source frontier_load_modules.sh


module reset
module load cpe/24.11
module load rocm/6.3.1
module load openblas/0.3.26
module load cray-hdf5
module load cmake
module load ninja
module load craype-accel-amd-gfx90a

export CC=mpicc
export CXX=mpicxx

export OPENBLAS_ROOT=${OLCF_OPENBLAS_ROOT}
export HDF5_ROOT=${OLCF_HDF5_ROOT}
export MAGMA_ROOT=/lustre/orion/world-shared/cph102/epd/opt/magma
export FFTW_PATH=${OLCF_FFTW_ROOT}
#export LD_PRELOAD=/opt/cray/pe/lib64/cce/libtcmalloc_minimal.so.1
export HDF5_ROOT=/lustre/orion/cph102/proj-shared/epd/spack/opt/spack/linux-sles15-zen3/rocmcc-6.0.0/hdf5-1.14.3-zb3ehzdumutzqowltqjhlzk6a3acywye
