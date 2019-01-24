# Initial cache list for Titan (Cray XK7)
#
# Usage: cmake /path/to/dca/source -C /path/to/this/file -D<option>=<value> ...

# Prevent CMake from searching for BLAS and LAPACK libraries.
# CC automatically links against them.
set(DCA_HAVE_LAPACK TRUE CACHE INTERNAL "")

# Use aprun for executing the tests.
set(TEST_RUNNER "aprun" CACHE STRING "Command for executing (MPI) programs.")

# Enable the GPU support.
option(DCA_WITH_CUDA "Enable GPU support." ON)

# Compile for Kepler compute architecture.
set(CUDA_GPU_ARCH "sm_35" CACHE STRING "Name of the *real* architecture to build for.")

# For the GPU support we also need MAGMA.
set(MAGMA_DIR $ENV{OLCF_MAGMA_ROOT} CACHE PATH
  "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")

# FFTW paths.
set(FFTW_INCLUDE_DIR $ENV{FFTW_INC} CACHE PATH "Path to fftw3.h.")
set(FFTW_LIBRARY $ENV{FFTW_DIR}/libfftw3.a CACHE FILEPATH "The FFTW3(-compatible) library.")
