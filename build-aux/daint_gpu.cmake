# Initial cache list for Piz Daint GPU (Cray XC50)
#
# Usage: cmake /path/to/dca/source -C /path/to/this/file -D<option>=<value> ...

# Prevent CMake from searching for BLAS and LAPACK libraries.
# CC automatically links against them.
set(DCA_HAVE_LAPACK TRUE CACHE INTERNAL "")

# Use srun for executing the tests.
set(TEST_RUNNER "srun" CACHE STRING "Command for executing (MPI) programs.")

# Enable the GPU support.
option(DCA_WITH_CUDA "Enable GPU support." ON)

# Compile for Tesla compute architecture.
set(CMAKE_CUDA_ARCHITECTURES 60 CACHE STRING "number of the CUDA compute architecture to build for.")  # default

#set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

set(FFTW_ROOT $ENV{FFTW_DIR}/.. CACHE PATH "Path to fftw3 library")

# For the GPU support we also need MAGMA.
# MAGMA has been installed with EasyBuild.
set(MAGMA_DIR $ENV{MAGMAROOT} CACHE PATH
  "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")

set(FFTW_INCLUDE_DIR $ENV{FFTW_ROOT}/include CACHE PATH "Path to fftw3.h.")
set(FFTW_LIBRARY $ENV{FFTW_ROOT}/lib/libfftw3.so CACHE FILEPATH "The FFTW3(-compatible) library.")

# Intel MKL flags
set(CMAKE_EXE_LINKER_FLAGS '-L/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl' CACHE INTERNAL "" FORCE)
