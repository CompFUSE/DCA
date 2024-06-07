################################################################################
# Author: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
#
# Checks for HIP and MAGMA and accordingly sets DCA_HAVE_HIP and DCA_HAVE_MAGMA.
# In addition, set DCA_GPU_LIBS.

# set(ROCM_ROOT
#   "/opt/rocm-4.2.0"
#     CACHE PATH "Root directory of ROCM")
# message(STATUS "ROCM_ROOT: ${ROCM_ROOT}")
# list(APPEND CMAKE_MODULE_PATH ${ROCM_ROOT}/hip/cmake)
# list(APPEND CMAKE_PREFIX_PATH ${ROCM_ROOT})
# find_package(HIP REQUIRED)
# find_package(hipblas REQUIRED)
# find_package(rocsolver REQUIRED)
# # architecture flags
# set(CMAKE_HIP_ARCHITECTURES
#   "gfx906,gfx908"
#   CACHE STRING "HIP architecture gfxXXX")
# list(APPEND HIP_HIPCC_FLAGS "-fPIC -ffast-math -O3")
# list(APPEND HIP_HIPCC_FLAGS "--amdgpu-target=${HIP_ARCH}")
# list(APPEND HIP_HIPCC_FLAGS "--gpu-max-threads-per-block=256")
# # warning suppression
# list(APPEND HIP_HIPCC_FLAGS "-Wno-vla")
# list(APPEND HIP_HIPCC_FLAGS "-Wno-deprecated-declarations")
# list(APPEND HIP_HIPCC_FLAGS "-Wno-unused-command-line-argument")
# list(APPEND HIP_HIPCC_FLAGS "-DHIP_PLATFORM_AMD")

# #-------------------------------------------------------------------
# #  set up ROCM compiler options and libraries
# #-------------------------------------------------------------------
if(DCA_WITH_HIP)
  if(${CMAKE_VERSION} VERSION_LESS "3.21.3")
    message(FATAL_ERROR "Compilation for HIP requires CMake 3.21.3 or later.")
  endif()
  set(ENABLE_HIP 1)
  message(STATUS "ROCM_ROOT: ${ROCM_ROOT}")

#-------------------------------------------------------------------
#  set up HIP compiler options
#-------------------------------------------------------------------
  set(CMAKE_MODULE_PATH "${ROCM_ROOT}/hip/cmake" "${ROCM_ROOT}/lib/cmake/hip" "${ROCM_ROOT}/lib/cmake/hipblas" "${ROCM_ROOT}/lib/cmake/rocthrust" ${CMAKE_MODULE_PATH})
  find_package(HIP REQUIRED)
  find_package(hipblas REQUIRED)
  find_package(hipsparse REQUIRED)
  find_package(rocsolver REQUIRED)
  find_package(rocthrust REQUIRED)

endif(DCA_WITH_HIP)

get_property(hipblas_include_dirs TARGET roc::hipblas PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
message("hipblas includes: ${hipblas_include_dirs}")


#set(CUDA_ARCHITECTURES "sm_60" CACHE STRING "Name of the real architecture to build for.")
set(MAGMA_ROOT "" CACHE PATH "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")

set(DCA_HAVE_HIP FALSE CACHE INTERNAL "")
set(DCA_HAVE_MAGMA FALSE CACHE INTERNAL "")
set(DCA_GPU_LIBS "" CACHE INTERNAL "")

include(CheckLanguage)
check_language(HIP)
if (CMAKE_HIP_COMPILER)
  enable_language(HIP)
  list(APPEND CMAKE_HIP_FLAGS "-fgpu-rdc")
  set(DCA_HAVE_HIP TRUE CACHE INTERNAL "")
  set(DCA_HAVE_GPU TRUE CACHE INTERNAL "")
  # Probably probably these should be public properties of the hip targets
  dca_add_haves_define(DCA_HAVE_HIP)
  dca_add_haves_define(DCA_HAVE_GPU)
  dca_add_haves_define(__HIP_PLATFORM_AMD__)
  list(APPEND DCA_GPU_LIBS hip::host roc::hipblas roc::hipsparse)
  set(DCA_HIP_PROPERTIES "CMAKE_HIP_ARCHITECTURES gfx908,gfx90a")
  set(CMAKE_HIP_STANDARD 17)
  list(APPEND HIP_HIPCC_FLAGS "-fPIC")
  list(APPEND HIP_HIPCC_FLAGS "-mno-unsafe-fp-atomics")
  list(APPEND HIP_HIPCC_FLAGS "-fgpu-default-stream=per-thread")
  list(APPEND HIP_HIPCC_FLAGS_DEBUG "--save-temps -g")

  # doesn't appear to work
  set(CMAKE_HIP_SOURCE_FILE_EXTENSIONS cu)
  message("Enabled HIP as a language")
  # NOTE: this is solved by dca_linking.cmake: dca_gpu_device_link()
  # alternative method (same issue)
  #file(GLOB_RECURSE CUDA_KERNELS_SRC ${PROJECT_SOURCE_DIR} *.cu)
  #set_source_files_properties(${CUDA_KERNELS_SRC} PROPERTIES LANGUAGE HIP)
# -ffast-math -O3")
# list(APPEND HIP_HIPCC_FLAGS "--amdgpu-target=${HIP_ARCH}")
# list(APPEND HIP_HIPCC_FLAGS "--gpu-max-threads-per-block=256")
# # warning suppression
# list(APPEND HIP_HIPCC_FLAGS "-Wno-vla")
# list(APPEND HIP_HIPCC_FLAGS "-Wno-deprecated-declarations")
# list(APPEND HIP_HIPCC_FLAGS "-Wno-unused-command-line-argument")
# list(APPEND HIP_HIPCC_FLAGS "-DHIP_PLATFORM_AMD")
endif()

#find_package(MAGMA
#  REQUIRED)

# # Find MAGMA.
find_package(MAGMA)
# find_library(MAGMA_LIBRARY
#   NAMES libmagma.so magma
#   HINTS ${MAGMA_DIR}/lib)
# find_path(MAGMA_INCLUDE_DIR magma.h HINTS ${MAGMA_DIR}/include)
# mark_as_advanced(MAGMA_LIBRARY MAGMA_INCLUDE_DIR)

if (MAGMA_LIBRARY AND MAGMA_INCLUDE_DIR)
  set(DCA_HAVE_MAGMA TRUE CACHE INTERNAL "")
  dca_add_haves_define(DCA_HAVE_MAGMA)
  # magma as of 2.2.0 is setup to build with openmp
  # if FindOpenMP.cmake finds it.
  # This can lead to link problems for us since
  # we don't otherwise have anything to do with openmp.
  # I have built magma without openmp for
  # CI. But if you naively use a random systems
  # magma expect to have a link error.
  target_link_libraries(magma::magma INTERFACE roc::hipblas roc::hipsparse LAPACK::LAPACK BLAS::BLAS)
  target_link_libraries(magma::sparse INTERFACE magma::magma)
  list(APPEND DCA_GPU_LIBS ${MAGMA_LIBRARY} roc::hipsparse roc::hipblas)
endif()
