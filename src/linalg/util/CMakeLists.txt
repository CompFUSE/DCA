# cuda and cublas utils

if(DCA_HAVE_GPU)
  set(srcs util_gpublas.cpp error_gpuBLAS.cpp info_gpu.cpp error_gpu.cpp)
  add_library(gpu_utils STATIC ${srcs})
  dca_gpu_runtime_link(gpu_utils)
  dca_gpu_blas_link(gpu_utils)
  target_link_libraries(gpu_utils PUBLIC magma::magma)
endif()

