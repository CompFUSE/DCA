#!/bin/bash
#
# Loads all modules that are required to build DCA++ on ORNL's Summit supercomputer.
# A reset is done at the beginning to restore to the default programming environment on Summit.
#
# Usage: source summit_load_modules.sh

module reset
module load gcc/9.3.0
module load cuda/11.4.0
module load magma/2.6.1
module load hdf5
module load fftw
module load cmake
module load netlib-lapack
module load essl

