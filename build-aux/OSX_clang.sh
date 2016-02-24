################################################################################
# Environment setup for OSX + clang
################################################################################

# Mock EasyBuild.

# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
LOCAL_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export EBROOTNFFT=$LOCAL_DIR/../libs
export EBROOTSPGLIB=$LOCAL_DIR/../libs
export EBROOTGTEST=$LOCAL_DIR/../libs/gmock-1.7.0/gtest

# Compilers
export CXX=/opt/local/bin/mpicxx-mpich-clang36
export CC=/opt/local/bin/mpicc-mpich-clang36
