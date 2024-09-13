# BLAS

if(DCA_HAVE_GPU)
  add_library(blas_kernels kernels_gpu.cu)
  target_link_libraries(blas_kernels PUBLIC magma::magma)
  dca_gpu_device_link(blas_kernels)
  dca_gpu_blas_link(blas_kernels)
  dca_gpu_runtime_link(blas_kernels)
endif()
