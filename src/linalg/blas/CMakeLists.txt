# BLAS

if(DCA_HAVE_GPU)
  add_library(blas_kernels kernels_gpu.cu)
  set_source_files_properties(blas_kernels kernels_gpu.cu PROPERTIES LANGUAGE HIP)
  dca_gpu_device_link(blas_kernels)
  dca_gpu_blas_link(blas_kernels)
  dca_gpu_runtime_link(blas_kernels)
  # set_target_properties( blas_kernels
  #                      PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  #                    target_compile_definitions(blas_kernels PRIVATE DCA_HAVE_CUDA)
  #                    target_link_libraries(blas_kernels PRIVATE CUDA::cublas)
endif()
