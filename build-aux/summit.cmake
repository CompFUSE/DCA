# Initial cache list for Summit
#
# Usage: cmake -C /path/to/this/file /path/to/DCA/source -D<option>=<value> ...

set(TEST_RUNNER "jsrun" CACHE STRING "Command for executing (MPI) programs.")
set(MPIEXEC_NUMPROC_FLAG "--np" CACHE STRING
  "Flag used by TEST_RUNNER to specify the number of processes.")
set(MPIEXEC_PREFLAGS "-g1 --smpiargs=none" CACHE STRING
  "Flags to pass to TEST_RUNNER directly before the executable to run.")

option(DCA_WITH_CUDA "Enable GPU support." ON)
set(CUDA_GPU_ARCH "sm_70" CACHE STRING "Name of the *real* architecture to build for.")
# Summit's static CUDA runtime is bugged.
option(CUDA_USE_STATIC_CUDA_RUNTIME OFF)

set(MAGMA_DIR $ENV{OLCF_MAGMA_ROOT} CACHE PATH
  "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")

set(FFTW_DIR $ENV{OLCF_FFTW_ROOT} CACHE PATH
  "Path to the FFTW3 installation directory. Hint for CMake to find FFTW3.")
