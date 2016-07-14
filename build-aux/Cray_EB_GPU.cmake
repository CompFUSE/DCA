# Initial cache list for Cray systems where the external libraries have been installed with
# EasyBuild
#
# GPU version
#
# Usage: cmake -C path/to/this/file/Cray_EB_GPU.cmake ...

option(DCA_WITH_CUDA "Enable CUDA support." ON)

# Include CPU cache list.
include(${CMAKE_CURRENT_LIST_DIR}/Cray_EB_CPU.cmake)

# For GPU support we also need MAGMA.
set(MAGMA_DIR $ENV{EBROOTMAGMA} CACHE FILEPATH "Path to MAGMA installation directory.")
mark_as_advanced(MAGMA_DIR)

# Compile for Kepler compute architecture.
set(CUDA_GPU_ARCH "sm_35" CACHE STRING "" FORCE)
