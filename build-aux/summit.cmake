# Initial cache list for Summit
#
# Usage: cmake -C /path/to/this/file /path/to/DCA/source -D<option>=<value> ...

# Prevent CMake from searching for BLAS and LAPACK libraries.
# IBM's ESSL (preferred) and LAPACK will be searched manually.
set(DCA_HAVE_LAPACK TRUE CACHE INTERNAL "")
set(LAPACK_LIBRARIES $ENV{OLCF_ESSL_ROOT}/lib64/libessl.so;$ENV{OLCF_NETLIB_LAPACK_ROOT}/lib64/liblapack.so CACHE FILEPATH " ")
message("LAPACK_LIBRARIES link:" ${LAPACK_LIBRARIES})

# Use jsrun for executing the tests.
# The flag "--smpiargs=none" is needed to execute tests with no MPI functionalities.
set(TEST_RUNNER "jsrun" CACHE STRING "Command for executing (MPI) programs.")
set(MPIEXEC_NUMPROC_FLAG "-a" CACHE STRING
  "Flag used by TEST_RUNNER to specify the number of processes.")
set(MPIEXEC_PREFLAGS "-n 1 -g 1 --smpiargs=none" CACHE STRING
  "Flags to pass to TEST_RUNNER directly before the executable to run.")

# Enable the GPU support.
option(DCA_WITH_CUDA "Enable GPU support." ON)

# Compile for Volta compute architecture.
set(CUDA_GPU_ARCH "sm_70" CACHE STRING "Name of the *real* architecture to build for.")

# Summit's static CUDA runtime is bugged.
option(CUDA_USE_STATIC_CUDA_RUNTIME OFF)

# For the GPU support we also need MAGMA.
set(MAGMA_DIR $ENV{OLCF_MAGMA_ROOT} CACHE PATH
  "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")

 # FFTW paths.
set(FFTW_DIR $ENV{OLCF_FFTW_ROOT} CACHE PATH
  "Path to the FFTW3 installation directory. Hint for CMake to find FFTW3.")
set(FFTW_INCLUDE_DIR $ENV{OLCF_FFTW_ROOT}/include CACHE PATH "Path to fftw3.h.")
set(FFTW_LIBRARY $ENV{OLCF_FFTW_ROOT}/lib/libfftw3.so CACHE FILEPATH "The FFTW3(-compatible) library.")
