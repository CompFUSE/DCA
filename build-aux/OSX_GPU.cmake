################################################################################
# Initial cache list for UNIX 
# Usage: cmake -C path/to/UNIX_clang.cmake ..
################################################################################

# call common option
include(${CMAKE_CURRENT_SOURCE_DIR}/../build-aux/UNIX_clang.cmake)

set(CUDA_GPU_ARCH "compute_50" CACHE STRING "gpu architecture" FORCE)
# Change dafault cache entries.
option(DCA_GPU_SUPPORT "Enable GPU support." ON)
option(DCA_PTHREADS    "Enable pthreads"     ON)
set(MAGMA_LIBRARY "/usr/local/lib/libmagma.a" CACHE FILEPATH "Location of magma lib" FORCE)
set(CUDA_HOST_COMPILER "/usr/bin/gcc-4.9"
  CACHE STRING "Host side compiler used by NVCC")
#set(CUDA_NVCC_FLAGS"${CUDA_NVCC_FLAGS}")
set(CUDA_SPARSE_LIBRARY "/usr/local/cuda-7.5/targets/x86_64-linux/lib/libcusparse.so"  CACHE FILEPATH "libcusparse path" FORCE)
mark_as_advanced(CUDA_GPU_ARCH CUDA_HOST_COMPILER) 
