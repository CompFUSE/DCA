if(DCA_HAVE_CUDA)
  CUDA_ADD_LIBRARY(blas_kernels kernels_gpu.cu)
endif()

