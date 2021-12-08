# BLAS

if(DCA_HAVE_GPU)
  add_library(blas_kernels kernels_gpu.cu)
  target_compile_definitions(blas_kernels PUBLIC DCA_HAVE_CUDA)
  dca_gpu_device_link(blas_kernels)
  dca_gpu_blas_link(blas_kernels)
  dca_gpu_runtime_link(blas_kernels)
endif()
