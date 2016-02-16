################################################################################
# Environment setup for OSX + clang
################################################################################

# Mock EasyBuild.

# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
LOCAL_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo ${LOCAL_DIR}
export EBROOTNFFT=$LOCAL_DIR/../libs
export EBROOTSPGLIB=$LOCAL_DIR/../libs
export EBROOTGTEST=$LOCAL_DIR/../libs/gmock-1.7.0/gtest

# MPI compilers
export MPI_CXX_COMPILER=/opt/local/bin/mpicxx-mpich-clang36
export MPI_C_COMPILER=/opt/local/bin/mpicc-mpich-clang36
