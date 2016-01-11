################################################################################
# Handles GPU support.
#
# Description.
################################################################################

option(DCA_GPU_SUPPORT            "Enable GPU support."        OFF)
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
    PATHS $ENV{EBROOTMAGMA}/lib
    NO_DEFAULT_PATH
    )

  set(DCA_GPU_LIBRARIES
    ${MAGMA_LIBRARY}
    ${CUDA_CUDA_LIBRARY}
    )
  # message("DCA_GPU_LIBRARIES: ${DCA_GPU_LIBRARIES}")
  
  CUDA_INCLUDE_DIRECTORIES(
    "${CMAKE_SOURCE_DIR}/src/comp_library/LIN_ALG"
    "${CMAKE_SOURCE_DIR}/src/phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker"
    "${CMAKE_SOURCE_DIR}/src/phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_N_matrix_routines"
    "${CMAKE_SOURCE_DIR}/src/phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_G_matrix_routines"
    "${CMAKE_SOURCE_DIR}/src/phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_G0_matrix_routines"
    "$ENV{EBROOTMAGMA}/include"
    )

else()
  # Source code replacements for CPU usage.
  set(DCA_CUDA_FUNCTION       "void print_device_info(){}")
  set(DCA_INITIALIZE_MAGMA_0  "")
  set(DCA_INITIALIZE_MAGMA_1  "")
  set(DCA_LIN_ALG_DEVICE_TYPE "LIN_ALG::CPU")

endif()
