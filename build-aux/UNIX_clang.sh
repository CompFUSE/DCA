#!/bin/bash
################################################################################
# Environment setup for UNIX.
# Remember to run as "source <this script>" before invoking cmake
################################################################################

# Mock EasyBuild.
export EBROOTNFFT=/usr/local
export EBROOTSPGLIB=/usr/local
export EBROOTGTEST=/home/giovanni/Programs/DCA/libs/gmock-1.7.0/gtest/

# Compilers
export CXX=/usr/lib/mpic++
export CC=/usr/lib/mpicc
