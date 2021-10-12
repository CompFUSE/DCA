# BLAS

if(DCA_HAVE_GPU)
  add_library(blas_kernels kernels_gpu.cu)
  if(DCA_HAVE_HIP)
    set_source_files_properties(blas_kernels kernels_gpu.cu PROPERTIES LANGUAGE HIP)
  endif()
  dca_gpu_device_link(blas_kernels)
  dca_gpu_blas_link(blas_kernels)
  dca_gpu_runtime_link(blas_kernels)
endif()
