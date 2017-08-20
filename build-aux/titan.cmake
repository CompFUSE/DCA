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

# Set the path to the CUDA toolkit since CMake cannot find it.
set(CUDA_TOOLKIT_ROOT_DIR "/opt/nvidia/cudatoolkit7.5/7.5.18-1.0502.10743.2.1"
  CACHE PATH "Path to the CUDA Toolkit." FORCE)

# For the GPU support we also need MAGMA.
# MAGMA has been installed with EasyBuild.
set(MAGMA_DIR $ENV{EBROOTMAGMA} CACHE PATH
  "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")
