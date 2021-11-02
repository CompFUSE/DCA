#!/bin/bash
#
# Loads all modules that are required to build DCA++ on ORNL's Spock cluster.
# A reset is done at the beginning to restore to the default programming environment on Spock.
# This is for development only at this point.
#
# Usage: source spock_load_cray_modules.sh

module reset
module load PrgEnv-cray
module load cmake/3.21
module load ninja
module load magma/2.6.1
module load fftw
module load cray-hdf5/1.12.0.6
module load cray-fftw/3.3.8.10
module load cray-libsci/21.06.1.1
module load cray-mpich/8.1.7
module load cray-python/3.8.5.1
module load rocm/4.2.0 # needed for magma 
