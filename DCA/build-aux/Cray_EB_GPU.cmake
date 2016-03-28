################################################################################
# Initial cache list for Cray systems where the external libraries are installed
# with EasyBuild.
# GPU version.
# Usage: cmake -C path/to/this/file/Cray_EB_GPU.cmake ...
################################################################################

# Include CPU cache list.
include($ENV{DCA_SOURCE}/build-aux/Cray_EB_CPU.cmake)

# For GPU support we also need MAGMA.
set(MAGMA_DIR $ENV{EBROOTMAGMA}
  CACHE FILEPATH "Path to MAGMA installation directory.")
mark_as_advanced(MAGMA_DIR)

# Workaround for CMake 3, otherwise FindCUDA sets CUDA_HOST_COMPILER
# incorrectly.
set(CUDA_HOST_COMPILER "${CMAKE_C_COMPILER}"
  CACHE FILEPATH "Host side compiler used by NVCC")  

# Compile for Kepler compute architecture.
set(CUDA_GPU_ARCH "sm_35" CACHE STRING "" FORCE)

# Change dafault cache entries.
option(DCA_GPU_SUPPORT "Enable GPU support." ON)
