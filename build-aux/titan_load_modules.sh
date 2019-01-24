#!/bin/bash
#
# Loads all modules that are required to build DCA++ on ORNL's Titan supercomputer.
# Assumes the default environment that one gets when logging into Titan.
#
# Usage: source titan_load_modules.sh

module unload PrgEnv-pgi
module load PrgEnv-gnu
module load cray-hdf5
module load cray-fftw
module load cudatoolkit
module load magma/2.4.0
module load cmake3
