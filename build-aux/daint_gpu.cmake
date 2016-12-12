# Initial cache list for Piz Daint XC50 (with GPU)
#
# Usage: cmake -C /path/to/this/file [other options] /path/to/source/dir

option(DCA_WITH_CUDA "Enable CUDA support." ON)

# Include the multicore cache list.
include(${CMAKE_CURRENT_LIST_DIR}/Cray_EB_CPU.cmake)

# For GPU support we also need MAGMA.
set(MAGMA_DIR $ENV{EBROOTMAGMA} CACHE PATH "Path to MAGMA installation directory.")

# Compile for Tesla compute architecture.
set(CUDA_GPU_ARCH "sm_60" CACHE STRING "" FORCE)

set(CUDA_TOOLKIT_ROOT_DIR "/opt/nvidia/cudatoolkit8.0/8.0.44_GA_2.2.7_g4a6c213-2.1"
  CACHE PATH "Path to the CUDA Toolkit." FORCE)
