###############################################################################
# Initial cache template for local machines
# Usage: cmake -C path/to/<this file> ..
################################################################################
include(local_CPU.cmake)

#magma library directory. The user is supposed to specify this.
set(DCA_LIBDIR_MAGMA "" CACHE FILEPATH "MAGMA directory")

#turn GPU support on
set(DCA_GPU_SUPPORT ON FORCE)

# Workaround for CMake 3, otherwise FindCUDA sets CUDA_HOST_COMPILER incorrectly.
#TODO try with new FindCuda
set(CUDA_HOST_COMPILER "${CMAKE_C_COMPILER}"
        CACHE FILEPATH "Host side compiler used by NVCC")
include(findCuda)

# Compile for Kepler compute architecture.
set(CUDA_GPU_ARCH "sm_35" CACHE STRING "" FORCE)
