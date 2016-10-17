if(DCA_HAVE_CUDA)
  CUDA_ADD_LIBRARY(lapack_kernels laset_gpu.cu multiply_diagonal_gpu.cu)
endif()

