################################################################################
# Environment setup for OSX + gcc
################################################################################

# Mock EasyBuild.
export EBROOTNFFT=$HOME/PhD/dca_ethz/libs
export EBROOTSPGLIB=$HOME/PhD/dca_ethz/libs
export EBROOTGTEST=$HOME/PhD/dca_ethz/libs/gmock-1.7.0/gtest

# Compilers
export CXX=/opt/local/bin/mpicxx-mpich-gcc5
export CC=/opt/local/bin/mpicc-mpich-gcc5
