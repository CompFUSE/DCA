################################################################################
# Initial cache list for OSX + gcc
# Usage: cmake -C path/to/OSX_gcc.cmake ...
################################################################################

# MPIEXEC stuff for executing parallel tests.
# Note: CXX and CC must be set as environment variables to the corresponding MPI
#       C++ and C compiler wrappers.
set(MPIEXEC "/opt/local/bin/mpiexec-mpich-gcc5"
  CACHE FILEPATH "Executable for running MPI programs.")
set(MPIEXEC_NUMPROC_FLAG "-np"
  CACHE FILEPATH
  "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes.")
