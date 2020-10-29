# Initial cache list for cades
#
# Building on this cluster is very brittle due to slurm and bad system level modules?
# Centos 7 in general?
#
# Spack generated hdf5 and magma seemed problematic so both are hand built.
#
# Don't expect this to work at all without sourcing
# build-aux/cades_load_modules.sh
#
# Usage: cmake -C /path/to/this/file /path/to/DCA/source -D<option>=<value> -GNinja ...

# Use srun for executing the tests.
set(TEST_RUNNER "srun" CACHE STRING "Command for executing (MPI) programs.")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING
  "Flag used by TEST_RUNNER to specify the number of processes.")
# Use 1 GPU and 64G memory per test process passing all tests will require running on 4 P100 nodes
# i.e.
# salloc -A ccsd -p gpu_p100 --nodes=4 --mem=180G --exclusive --gres=gpu:2 -t 00:30:00
set(MPIEXEC_PREFLAGS "--mem=64G --gpus-per-task=1" CACHE STRING
  "Flags to pass to TEST_RUNNER directly before the executable to run.")

# these aren't needed on cades.
set(SMPIARGS_FLAG_NOMPI "" CACHE STRING
  "Spectrum MPI argument list flag for serial tests.")
# Let's keep this option in case we need it again in the future.
set(SMPIARGS_FLAG_MPI "" CACHE STRING "Spectrum MPI argument list flag for MPI tests.")

# Enable the GPU support.
option(DCA_WITH_CUDA "Enable GPU support." ON)
option(DCA_WITH_CUDA_AWARE_MPI "Enable CUDA aware MPI." OFF)

set(CUDA_TOOLKIT_ROOT_DIR $ENV{CUDA_DIR} CACHE PATH "path to CUDA toolkit")

# Compile for Volta compute architecture.
set(CUDA_GPU_ARCH "sm_60" CACHE STRING "Name of the *real* architecture to build for.")

# Summit's static CUDA runtime is bugged.
option(CUDA_USE_STATIC_CUDA_RUNTIME OFF)

# For the GPU support we also need MAGMA.
set(MAGMA_DIR $ENV{MAGMA_DIR} CACHE PATH
  "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")

# FFTW paths.
set(FFTW_INCLUDE_DIR $ENV{FFTW_DIR}/include CACHE PATH "Path to fftw3.h.")
set(FFTW_LIBRARY $ENV{FFTW_DIR}/lib/libfftw3.so CACHE FILEPATH "The FFTW3(-compatible) library.")

# HDF5 paths
set(HDF5_ROOT $ENV{HD5_DIR})
set(HDF5_INCLUDE_DIRS $ENV{HDF5_DIR}/include CACHE PATH "Path to hdf5 includes")
set(HDF5_LIBRARIES "$ENV{HDF5_DIR}/lib/libhdf5_cpp.a;$ENV{HDF5_DIR}/lib/libhdf5.a" CACHE FILEPATH "The hdf5 libraries")

option(DCA_WITH_TESTS_FAST "Fast minimal tests" ON)

#required by dependencies but not picked up by cmake for whatever reason.
set(CMAKE_EXE_LINKER_FLAGS "-ldl -fopenmp" CACHE STRING "additional linking arguments needed")
