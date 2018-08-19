# BLAS

if(DCA_HAVE_CUDA)
  CUDA_ADD_LIBRARY(blas_kernels kernels_gpu.cu)
  target_compile_definitions(blas_kernels PRIVATE DCA_HAVE_CUDA)
endif()
