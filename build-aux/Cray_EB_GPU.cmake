# Initial cache list for Cray systems where the external libraries have been installed with
# EasyBuild
#
# GPU version
#
# Usage: cmake -C /path/to/this/file [other options] /path/to/source/dir

option(DCA_WITH_CUDA "Enable CUDA support." ON)

# Include the multicore cache list.
include(${CMAKE_CURRENT_LIST_DIR}/Cray_EB_CPU.cmake)

# For GPU support we also need MAGMA.
set(MAGMA_DIR $ENV{EBROOTMAGMA} CACHE PATH "Path to MAGMA installation directory.")

# Compile for Kepler compute architecture.
# TODO: Need sm_60 for new Piz Daint!
set(CUDA_GPU_ARCH "sm_35" CACHE STRING "" FORCE)
