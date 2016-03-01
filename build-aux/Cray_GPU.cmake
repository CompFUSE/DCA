################################################################################
# Initial cache list for Cray systems with GPU (Piz Daint, Titan)
# Usage: cmake -C path/to/Cray_GPU.cmake ...
################################################################################


#include common definitioons from CPU setup
include(${CMAKE_CURRENT_SOURCE_DIR}/Cray_CPU.cmake)
# The C++ compile wrapper CC already includes and links to these libraries.
# No need to look for them.
set(DCA_HDF5_IMPLICIT   TRUE CACHE INTERNAL "")
set(DCA_LAPACK_IMPLICIT TRUE CACHE INTERNAL "")
set(DCA_FFTW_IMPLICIT   TRUE CACHE INTERNAL "")

# Workaround for CMake 3, otherwise FindCUDA sets CUDA_HOST_COMPILER incorrectly.
set(CUDA_HOST_COMPILER "${CMAKE_C_COMPILER}"
  CACHE FILEPATH "Host side compiler used by NVCC")  

# Compile for Kepler compute architecture.
set(CUDA_GPU_ARCH "sm_35" CACHE STRING "" FORCE)

# Change dafault cache entries.
set(DCA_GPU_SUPPORT ON FORCE)

