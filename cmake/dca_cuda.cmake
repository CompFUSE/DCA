################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Checks for CUDA and MAGMA and accordingly sets DCA_HAVE_CUDA and DCA_HAVE_MAGMA.
# In addition, set DCA_CUDA_LIBS.

set(CUDA_GPU_ARCH "sm_60" CACHE STRING "Name of the real architecture to build for.")
set(MAGMA_DIR "" CACHE PATH "Path to the MAGMA installation directory. Hint for CMake to find MAGMA.")

set(DCA_HAVE_CUDA FALSE CACHE INTERNAL "")
set(DCA_HAVE_MAGMA FALSE CACHE INTERNAL "")
set(DCA_CUDA_LIBS "" CACHE INTERNAL "")

set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)

# Find CUDA.
find_package(CUDA REQUIRED)

if (CUDA_FOUND)
  # set(DCA_HAVE_CUDA TRUE CACHE INTERNAL "")
  # dca_add_haves_define(DCA_HAVE_CUDA)
  list(APPEND DCA_CUDA_LIBS ${CUDA_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY})
  CUDA_INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})
  set(CUDA_SEPARABLE_COMPILATION ON)
endif()

# Find MAGMA.
find_library(MAGMA_LIBRARY
  NAMES libmagma.a magma
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
  list(APPEND DCA_CUDA_LIBS ${MAGMA_LIBRARY} ${CUDA_cusparse_LIBRARY})
  CUDA_INCLUDE_DIRECTORIES(${MAGMA_INCLUDE_DIR})
endif()

# At the moment the GPU code requires MAGMA. Therefore we set DCA_HAVE_CUDA to true, only if both
# CUDA and MAGMA have been found.
if (CUDA_FOUND AND DCA_HAVE_MAGMA)
  set(DCA_HAVE_CUDA TRUE CACHE INTERNAL "")
  dca_add_haves_define(DCA_HAVE_CUDA)
endif()
