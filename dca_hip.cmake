################################################################################
# Author: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
#
# Checks for HIP and MAGMA and accordingly sets DCA_HAVE_HIP and DCA_HAVE_MAGMA.
# In addition, set DCA_HIP_LIBS.

set(ROCM_ROOT
    "/opt/rocm"
    CACHE PATH "Root directory of ROCM")
message(STATUS "ROCM_ROOT: ${ROCM_ROOT}")
list(APPEND CMAKE_MODULE_PATH ${ROCM_ROOT}/hip/cmake)
list(APPEND CMAKE_PREFIX_PATH ${ROCM_ROOT})
find_package(HIP REQUIRED)
find_package(hipblas REQUIRED)
find_package(rocsolver REQUIRED)
# architecture flags
set(HIP_ARCH
  "gfx906,gfx908"
  CACHE STRING "HIP architecture gfxXXX")
list(APPEND HIP_HIPCC_FLAGS "-fPIC -ffast-math -O3")
list(APPEND HIP_HIPCC_FLAGS "--amdgpu-target=${HIP_ARCH}")
list(APPEND HIP_HIPCC_FLAGS "--gpu-max-threads-per-block=256")
# warning suppression
list(APPEND HIP_HIPCC_FLAGS "-Wno-vla")
list(APPEND HIP_HIPCC_FLAGS "-Wno-deprecated-declarations")
list(APPEND HIP_HIPCC_FLAGS "-Wno-unused-command-line-argument")

#-------------------------------------------------------------------
#  set up ROCM compiler options and libraries
#-------------------------------------------------------------------
if(ENABLE_ROCM)
  message(STATUS "ROCM_ROOT: ${ROCM_ROOT}")
  add_library(ROCM::libraries INTERFACE IMPORTED)
  # temporarily put rocsolver rocrand here for convenience, should be moved to Platforms.
  set_target_properties(ROCM::libraries PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${ROCM_ROOT}/include"
                                                   INTERFACE_LINK_LIBRARIES "-L${ROCM_ROOT}/lib;-lrocsolver;-lrocrand")
endif(ENABLE_ROCM)

#-------------------------------------------------------------------
#  set up HIP compiler options
#-------------------------------------------------------------------
if(ENABLE_HIP)
  if(NOT ENABLE_ROCM)
    message(FATAL_ERROR "ROCM is required to use HIP. Please set ENABLE_ROCM=ON.")
  endif()
  set(CMAKE_MODULE_PATH "${ROCM_ROOT}/hip/cmake" ${CMAKE_MODULE_PATH})
  find_package(HIP REQUIRED)

  add_library(HIP::HIP INTERFACE IMPORTED)
  # temporarily put hipsparse hipblas here for convenience, should be moved to Platforms.
  set_target_properties(
    HIP::HIP PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${ROCM_ROOT}/include" INTERFACE_COMPILE_DEFINITIONS "ENABLE_HIP"
                        INTERFACE_LINK_LIBRARIES "-L${ROCM_ROOT}/lib;-lhipsparse;-lhipblas;-lamdhip64")
endif(ENABLE_HIP)

#set(CUDA_ARCHITECTURES "sm_60" CACHE STRING "Name of the real architecture to build for.")
set(MAGMA_DIR "" CACHE PATH "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")

set(DCA_HAVE_GPU TRUE CACHE INTERNAL "")
set(DCA_HAVE_HIP FALSE CACHE INTERNAL "")
set(DCA_HAVE_MAGMA FALSE CACHE INTERNAL "")
set(DCA_HIP_LIBS "" CACHE INTERNAL "")

include(CheckLanguage)
check_language(HIP)
if (CMAKE_HIP_COMPILER)
  enable_language(HIP)
  set(DCA_HAVE_HIP TRUE CACHE INTERNAL "")
  set(DCA_HAVE_GPU TRUE CACHE INTERNAL "")
  dca_add_haves_define(DCA_HAVE_HIP)
  list(APPEND DCA_HIP_LIBS
    HIP::HIP ROCM::libraries)
  set(DCA_HIP_PROPERTIES "CMAKE_HIP_ARCHITECTURES gfx906,gfx908")
  set(CMAKE_HIP_STANDARD 17)
endif()

#find_package(MAGMA
#  REQUIRED)

# # Find MAGMA.
find_library(MAGMA_LIBRARY
  NAMES libmagma.so magma
  HINTS ${MAGMA_DIR}/lib)
find_path(MAGMA_INCLUDE_DIR magma.h HINTS ${MAGMA_DIR}/include)
mark_as_advanced(MAGMA_LIBRARY MAGMA_INCLUDE_DIR)

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
endif()

# At the moment the GPU code requires MAGMA. Therefore we set DCA_HAVE_CUDA to true, only if both
# CUDA and MAGMA have been found.
if (HIP_FOUND AND DCA_HAVE_MAGMA)
  set(DCA_HAVE_HIP TRUE CACHE INTERNAL "")
  dca_add_haves_define(DCA_HAVE_HIP)
endif()
