#!/bin/bash
#
# Loads all modules that are required to build DCA++ on ORNL's Summit supercomputer.
# A reset is done at the beginning to restore to the default programming environment on Summit.
#
# Usage: source adios2_summit_load_modules.sh

module reset
module load gcc/8.1.1
module load cuda/11.1.1
module load hdf5
module load fftw
module load cmake
module load netlib-lapack
module load essl


export MAGMA_DIR=/ccs/home/pdoak/opt/magma_cuda11
export CC=mpicc
export CXX=mpicxx
module load bzip2/1.0.6

# Right now this assumes you have a build of ADIOS2 that
# must be passed via ADIOS2_DIR.  Eventually a working module will be required.
#
# cmake -C ../build-aux/summit.cmake -DDCA_WITH_MPI=ON -DDCA_WITH_TESTS_FAST=ON -DADIOS2_DIR=/gpfs/alpine/proj-shared/cph102/epd/ADIOS2_summit -DCMAKE_BUILD_TYPE=Debug -GNinja ..
