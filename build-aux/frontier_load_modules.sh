#!/bin/bash
#
# Loads all modules that are required to build DCA++ on ORNL's Frontier.
# A reset is done at the beginning to restore to the default programming environment on Frontier.
# This is for development only at this point.
#
# Usage: source frontier_load_modules.sh

module reset
module load gcc/12.2.0
module load rocm/5.7.0
module load ninja

# After 2 weeks of digging through opaque linking and runtime errors,
# I have concluded that cray-libsci causes such a mess
# that it's much easier to compile your own openblas
# and magma rather than fuss with it.  I did the latter in 1 day.
module unload cray-libsci

export CC=mpicc
export CXX=mpicxx
