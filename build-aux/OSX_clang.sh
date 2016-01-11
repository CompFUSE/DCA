################################################################################
# Environment setup for OSX + clang
################################################################################

# Mock EasyBuild.
export EBROOTNFFT=$HOME/PhD/dca_ethz/libs
export EBROOTSPGLIB=$HOME/PhD/dca_ethz/libs
export EBROOTGTEST=$HOME/PhD/dca_ethz/libs/gmock-1.7.0/gtest

# Compilers
export CXX=/opt/local/bin/mpicxx-mpich-clang36
export CC=/opt/local/bin/mpicc-mpich-clang36
