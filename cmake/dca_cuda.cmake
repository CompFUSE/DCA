################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Checks for CUDA and MAGMA and accordingly sets DCA_HAVE_CUDA and DCA_HAVE_MAGMA.
#
# TODO: Check this file.

set(DCA_HAVE_CUDA FALSE CACHE INTERNAL "")
set(DCA_HAVE_MAGMA FALSE CACHE INTERNAL "")
set(DCA_CUDA_LIBS "" CACHE INTERNAL "")

# Path to libcuda.so on Daint and Titan. Needed by FindCUDA.
# TODO: Can we move this to the init-cache-file for daint/titan?
set(ENV{CUDA_LIB_PATH} "/opt/cray/nvidia/default/lib64")

# Find CUDA.
find_package(CUDA)

if (CUDA_FOUND)
  # set(DCA_HAVE_CUDA TRUE CACHE INTERNAL "")
  # dca_add_haves_define(DCA_HAVE_CUDA)

  set(CUDA_PROPAGATE_HOST_FLAGS OFF
    CACHE BOOL "Propage C/CXX_FLAGS and friends to the host compiler via -Xcompile")

  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=${CUDA_GPU_ARCH} -DNDEBUG")
  list(APPEND DCA_CUDA_LIBS ${CUDA_CUDA_LIBRARY} ${CUDA_cusparse_LIBRARY})
  CUDA_INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
endif()

# Find MAGMA.
find_library(MAGMA_LIBRARY magma HINTS ${MAGMA_DIR}/lib)
find_path(MAGMA_INCLUDE_DIR magma.h HINTS ${MAGMA_DIR}/include)
mark_as_advanced(MAGMA_LIBRARY MAGMA_INCLUDE_DIR)

if (MAGMA_LIBRARY AND MAGMA_INCLUDE_DIR)
  set(DCA_HAVE_MAGMA TRUE CACHE INTERNAL "")
  dca_add_haves_define(DCA_HAVE_MAGMA)

  # INTERNAL: When MAGMA is not required anymore for all the GPU code, we could remove/modify the
  #           next two lines.
  list(APPEND DCA_CUDA_LIBS ${MAGMA_LIBRARY})
  CUDA_INCLUDE_DIRECTORIES(${MAGMA_INCLUDE_DIR})
endif()

# At the moment the GPU code requires MAGMA. Therefore we set DCA_HAVE_CUDA to true, only if both
# CUDA and MAGMA have been found.
if (CUDA_FOUND AND DCA_HAVE_MAGMA)
  set(DCA_HAVE_CUDA TRUE CACHE INTERNAL "")
  dca_add_haves_define(DCA_HAVE_CUDA)
endif()
