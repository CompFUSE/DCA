#!/bin/bash
#
# Loads all modules that are required to build DCA++ on ORNL's Summit supercomputer.
# A reset is done at the beginning to restore to the default programming environment on Summit.
#
# Usage: source summit_load_modules.sh

module reset
module load gcc/10.2.0
module load cuda/11.1.1 # ldd shows magma is built with this cuda
module load magma/2.6.1
module load hdf5
module load fftw
module load cmake/3.20.2 # at least 3.20 is required
module load netlib-lapack
module load essl

