################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Checks for CUDA and MAGMA and accordingly sets DCA_HAVE_CUDA and DCA_HAVE_MAGMA.
# In addition, set DCA_GPU_LIBS.

set(CUDA_ARCHITECTURES "sm_70" CACHE STRING "Name of the real architecture to build for.")
set(MAGMA_DIR "" CACHE PATH "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")

set(DCA_HAVE_CUDA FALSE CACHE INTERNAL "")
set(DCA_HAVE_MAGMA FALSE CACHE INTERNAL "")
set(DCA_GPU_LIBS "" CACHE INTERNAL "")

set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)

# Find CUDA.
#find_package(CUDA REQUIRED)
include(CheckLanguage)

find_package(CUDAToolkit REQUIRED)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
  enable_language(CUDA)
  set(DCA_HAVE_CUDA TRUE CACHE INTERNAL "")
  set(DCA_HAVE_GPU TRUE CACHE INTERNAL "")
  dca_add_haves_define(DCA_HAVE_CUDA)
  dca_add_haves_define(DCA_HAVE_GPU)

  list(APPEND DCA_GPU_LIBS CUDA::cudart CUDA::cublas)
  set(DCA_CUDA_PROPERTIES "CMAKE_CUDA_ARCHITECTURES 70")
  list(APPEND CUDAFLAGS "--expt-relaxed-constexpr" ${DCA_CUDA_OPTIONS})
  set(CMAKE_CUDA_STANDARD 14)
  set(CVD_LAUNCHER "" CACHE INTERNAL "launch script for setting the Cuda visible devices.")
  # Use the following script for systems with multiple gpus visible from a rank.
  # set(CVD_LAUNCHER "test/cvd_launcher.sh" CACHE INTERNAL "")
endif()

#find_package(MAGMA
#  REQUIRED)

# # Find MAGMA.
find_package(MAGMA)

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
  target_link_libraries(magma::sparse INTERFACE magma::magma CUDA::cublas CUDA::cusparse)
  list(APPEND DCA_GPU_LIBS ${MAGMA_LIBRARY} CUDA::cusparse)
endif()
