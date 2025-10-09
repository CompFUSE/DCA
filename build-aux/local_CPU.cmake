# Template for initial cache list - CPU version (no GPU support)
#
# Modify for your local machine.
#
# Usage: cmake -C /path/to/this/file/local_CPU.cmake ...

# Location of external libraries
#
# SPRNG
set(SPRNG_DIR "/Projects/ETH_Reps/SPRNG/sprng5/"
  CACHE PATH "Path to SPRNG installation directory.")

# Google Test: This uses the provided version of Google Test in libs.
set(gtest_DIR "${CMAKE_CURRENT_LIST_DIR}/../libs/gmock-1.7.0/gtest/"
  CACHE PATH "Path to Google Test.")

# HDF5
# CMake 3.6 does not properly find HDF5 if HDF5_ROOT is not set.
set(HDF5_ROOT "/opt/homebrew/Cellar/hdf5/1.14.4.3/" CACHE PATH "Path to HDF5 installation directory.")
# set(HDF5_ROOT "/usr/local/anaconda3/pkgs/hdf5-1.10.6-hdbbcd12_0/" CACHE PATH "Path to HDF5 installation directory.")

# FFTW
# set(FFTW_INCLUDE_DIR "/usr/local/Cellar/fftw/3.3.10/include" CACHE PATH "Path to fftw3.h.")
# set(FFTW_LIBRARY "/usr/local/Cellar/fftw/3.3.10/lib/libfftw3.a" CACHE FILEPATH "Path to FFTW3 library.")
set(FFTW_INCLUDE_DIR "/opt/homebrew/Cellar/fftw/3.3.10_1/include" CACHE PATH "Path to fftw3.h.")
set(FFTW_LIBRARY "/opt/homebrew/Cellar/fftw/3.3.10_1/lib/libfftw3.a" CACHE FILEPATH "Path to FFTW3 library.")

mark_as_advanced(gtest_DIR SPRNG_DIR HDF5_ROOT)

# TEST_RUNNER must be set to the command to run an application.
# In particular, if CXX and CC are MPI C++ and C compiler wrappers, it is the command that runs an
# MPI application.
set(TEST_RUNNER "/usr/local/bin/mpiexec"
  CACHE FILEPATH "Command to run an application.")

set(MPIEXEC_NUMPROC_FLAG "-np"
  CACHE STRING "Flag used by TEST_RUNNER to specify the number of processes.")

mark_as_advanced(TEST_RUNNER MPIEXEC_NUMPROC_FLAG)

