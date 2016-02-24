################################################################################
# Initial cache list for OSX + clang
# Usage: cmake -C path/to/OSX_clang.cmake ...
################################################################################

# MPIEXEC stuff for executing parallel tests.
# Note: CXX and CC must be set as environment variables to the corresponding MPI
#       C++ and C compiler wrappers.
set(MPIEXEC "/opt/local/bin/mpiexec-mpich-clang36"
  CACHE FILEPATH "Executable for running MPI programs.")
set(MPIEXEC_NUMPROC_FLAG "-np"
  CACHE FILEPATH
  "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes.")
