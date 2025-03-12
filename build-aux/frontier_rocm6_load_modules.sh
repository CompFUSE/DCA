#!/bin/bash
#
# Loads all modules that are required to build DCA++ on ORNL's Frontier.
# A reset is done at the beginning to restore to the default programming environment on Frontier.
# This is for development only at this point.
#
# Usage: source frontier_load_modules.sh


module reset
module load rocm/6.0.0
module load ninja
module load cmake

export CC=mpicc
export CXX=mpicxx
export HDF5_ROOT=/lustre/orion/cph102/proj-shared/epd/spack/opt/spack/linux-sles15-zen3/rocmcc-6.0.0/hdf5-1.12.1-ajskwiaabdvgc36ozb6hzqnrwu2becha
export OPENBLAS_ROOT=/lustre/orion/cph102/proj-shared/epd/spack/opt/spack/linux-sles15-zen3/gcc-11.2.0/openblas-0.3.25-scaywvuh5zsm5u7smg54plj2oyf7nekv
export MAGMA_ROOT=/lustre/orion/cph102/proj-shared/epd/spack/opt/spack/linux-sles15-zen3/rocmcc-6.0.0/magma-master-rizw3ajkhfcq5cjutoykgkkv5hexftoz
export FFTW_PATH=/lustre/orion/cph102/proj-shared/epd/spack/opt/spack/linux-sles15-zen3/rocmcc-6.0.0/fftw-3.3.10-2mykijticsr5rfbyunax4zrwhhzcb7qm
#export LD_PRELOAD=/opt/cray/pe/lib64/cce/libtcmalloc_minimal.so.1
