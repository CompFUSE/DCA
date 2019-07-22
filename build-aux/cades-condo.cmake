# Initial cache list for CADES-CONDO
#
# Usage: cmake -C /path/to/this/file /path/to/DCA/source -D<option>=<value> ...

set(CMAKE_CXX_FLAGS "-fPIC -std=c++14" CACHE STRING "compiler flags")

set(FFTW_INCLUDE_DIR $ENV{FFTW_DIR}/include CACHE PATH "FFTW include path")
set(FFTW_LIBRARY $ENV{FFTW_DIR}/lib/libfftw3.a CACHE FILEPATH "FFTW libary")

set(TEST_RUNNER "mpirun" CACHE STRING "Command for executing (MPI) programs.")
set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING
  "Flag used by TEST_RUNNER to specify the number of processes.")
set(MPIEXEC_PREFLAGS "" CACHE STRING
  "Flags to pass to TEST_RUNNER directly before the executable to run.")

# Enable the GPU support.
option(DCA_WITH_CUDA "Enable GPU support." ON)
option(DCA_WITH_MPI "Enable MPI support." ON)
# Compile for Volta compute architecture.
set(CUDA_GPU_ARCH "sm_60" CACHE STRING "Name of the *real* architecture to build for.")

# For the GPU support we also need MAGMA.
#set(MAGMA_DIR $ENV{MAGMA_DIR} CACHE PATH)
#  "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")

set(CUDA_TOOLKIT_ROOT_DIR $ENV{CUDA_DIR})

option(DCA_WITH_TESTS_FAST "Enable fast tests" ON)

