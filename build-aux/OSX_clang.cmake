################################################################################
# Initial cache list for OSX + clang
# Usage: cmake -C path/to/OSX_clang.cmake ...
################################################################################
set(EBROOTNFFT "${CMAKE_SOURCE_DIR}/../libs" CACHE FILEPATH "NFFT location")
set(EBROOTSPGLIB  "${CMAKE_SOURCE_DIR}/../libs" CACHE FILEPATH "idem as before")
set(EBROOTGTEST "${CMAKE_CURRENT_SOURCE_DIR}/../libs/gmock-1.7.0/gtest/" CACHE FILEPATH "folder with gtest")
# Compilers for MPI
set(MPI_CXX_COMPILER /opt/local/bin/mpicxx-mpich-clang36)
set(MPI_C_COMPILER /opt/local/bin/mpicc-mpich-clang36)
#options and flags
option(DCA_GPU_SUPPORT "Disable GPU support." OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS};-std=c++11;-fopenmp")


# MPIEXEC stuff for executing parallel tests.
# Note: CXX and CC must be set as environment variables to the corresponding MPI
#       C++ and C compiler wrappers.
set(MPIEXEC "/opt/local/bin/mpiexec-mpich-clang36"
  CACHE FILEPATH "Executable for running MPI programs.")
set(MPIEXEC_NUMPROC_FLAG "-np"
  CACHE FILEPATH
  "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes.")
