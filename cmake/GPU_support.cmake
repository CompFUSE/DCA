################################################################################
# Handles GPU support.
#
# Description.
################################################################################

option(DCA_GPU_SUPPORT "Enable GPU support." OFF)
option(DCA_GPU_PINNED_HOST_MEMORY "Enable pinned host memory." OFF)
mark_as_advanced(DCA_GPU_PINNED_HOST_MEMORY)

if (DCA_GPU_SUPPORT)
  # Path to libcuda.so in Daint and Titan. Needed by FindCUDA.
  set(ENV{CUDA_LIB_PATH} "/opt/cray/nvidia/default/lib64")

  # FIXME: Only necessary when CXX_FLAGS include C++11.
  set(CUDA_PROPAGATE_HOST_FLAGS OFF
    CACHE BOOL "Propage C/CXX_FLAGS and friends to the host compiler via -Xcompile")
  
  find_package(CUDA REQUIRED)

  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=${CUDA_GPU_ARCH} -DNDEBUG")
  # message("${CUDA_NVCC_FLAGS}")

  add_definitions(-DUSE_GPU)

  if (DCA_GPU_PINNED_HOST_MEMORY)
    add_definitions(-DENABLE_PINNED_MEMORY_ALLOCATION)
  endif()
  
  # Source code replacements for GPU usage.
  set(DCA_CUDA_FUNCTION       "void print_device_info();")
  set(DCA_INITIALIZE_MAGMA_0  "void initialize_magma();")
  set(DCA_INITIALIZE_MAGMA_1  "initialize_magma();")
  set(DCA_LIN_ALG_DEVICE_TYPE "LIN_ALG::GPU")

  find_library(MAGMA_LIBRARY
    NAMES libmagma.a magma
    PATHS ${MAGMA_DIR}/lib
    NO_DEFAULT_PATH
    )
  mark_as_advanced(MAGMA_LIBRARY)

  set(DCA_GPU_LIBRARIES
    ${MAGMA_LIBRARY}
    ${CUDA_CUDA_LIBRARY}
    ${CUDA_cusparse_LIBRARY}
    )
  
  CUDA_INCLUDE_DIRECTORIES(
    "${PROJECT_SOURCE_DIR}/src"
    "${MAGMA_DIR}/include"
    )

else()
  # Source code replacements for CPU usage.
  set(DCA_CUDA_FUNCTION       "void print_device_info(){}")
  set(DCA_INITIALIZE_MAGMA_0  "")
  set(DCA_INITIALIZE_MAGMA_1  "")
  set(DCA_LIN_ALG_DEVICE_TYPE "LIN_ALG::CPU")

endif()
