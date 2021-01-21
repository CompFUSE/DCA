# BLAS

if(DCA_HAVE_CUDA)
  add_library(blas_kernels kernels_gpu.cu)
  set_target_properties( blas_kernels
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
                     target_compile_definitions(blas_kernels PRIVATE DCA_HAVE_CUDA)
                     target_link_libraries(blas_kernels PRIVATE CUDA::cublas)
endif()
