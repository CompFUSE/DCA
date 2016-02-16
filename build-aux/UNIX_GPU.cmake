################################################################################
# Initial cache list for UNIX 
# Usage: cmake -C path/to/UNIX_clang.cmake ..
################################################################################

# MPIEXEC stuff for executing parallel tests.
# Note: CXX and CC must be set as environment variables to the corresponding MPI
#       C++ and C compiler wrappers.

# Compile for Maxwell compute architecture.
set(CUDA_GPU_ARCH "compute_50" CACHE STRING "gpu architecture" FORCE)
# Change dafault cache entries.
option(DCA_GPU_SUPPORT "Enable GPU support." ON)
option(DCA_PTHREADS    "Enable pthreads"     ON)
option( MAGMA_LIBRARY "Locate magma lib" "/usr/local/lib/libmagma.so")
set(CUDA_HOST_COMPILER "/usr/bin/gcc-4.9"
  CACHE STRING "Host side compiler used by NVCC")
mark_as_advanced(CUDA_GPU_ARCH CUDA_HOST_COMPILER)

set(MPIEXEC "/usr/bin/mpirun"
  CACHE FILEPATH "Executable for running MPI programs.")
set(MPIEXEC_NUMPROC_FLAG "-np"
  CACHE FILEPATH
  "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes.")
