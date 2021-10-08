################################################################################
# Author: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
#
# Checks for HIP and MAGMA and accordingly sets DCA_HAVE_HIP and DCA_HAVE_MAGMA.
# In addition, set DCA_HIP_LIBS.

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
  set(ENABLE_HIP 1)
  message(STATUS "ROCM_ROOT: ${ROCM_ROOT}")

#-------------------------------------------------------------------
#  set up HIP compiler options
#-------------------------------------------------------------------
  set(CMAKE_MODULE_PATH "${ROCM_ROOT}/hip/cmake" ${CMAKE_MODULE_PATH})
  find_package(HIP REQUIRED)
  find_package(hipblas REQUIRED)
  find_package(hipsparse REQUIRED)
  find_package(rocsolver REQUIRED)

endif(DCA_WITH_HIP)

#set(CUDA_ARCHITECTURES "sm_60" CACHE STRING "Name of the real architecture to build for.")
set(MAGMA_ROOT "" CACHE PATH "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")

set(DCA_HAVE_HIP FALSE CACHE INTERNAL "")
set(DCA_HAVE_MAGMA FALSE CACHE INTERNAL "")
set(DCA_HIP_LIBS "" CACHE INTERNAL "")

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
  list(APPEND DCA_HIP_LIBS
    HIP::HIP ROCM::libraries)
  set(DCA_HIP_PROPERTIES "CMAKE_HIP_ARCHITECTURES gfx906,gfx908")
  set(CMAKE_HIP_STANDARD 17)
  list(APPEND HIP_HIPCC_FLAGS "-fPIC")
  # doesn't appear to work
  set(CMAKE_HIP_SOURCE_FILE_EXTENSIONS cu)
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
  list(APPEND DCA_HIP_LIBS ${MAGMA_LIBRARY} HIP::sparse)
  target_link_libraries(magma::magma INTERFACE LAPACK::LAPACK BLAS::BLAS)
endif()

# At the moment the GPU code requires MAGMA. Therefore we set DCA_HAVE_HIP to true, only if both
# HIP and MAGMA have been found.
if (HIP_FOUND AND DCA_HAVE_MAGMA)
  set(DCA_HAVE_HIP TRUE CACHE INTERNAL "")
  MESSAGE("DCA_HAVE_HIP set with dca_add_haves_define")
  dca_add_haves_define(DCA_HAVE_HIP)
endif()
