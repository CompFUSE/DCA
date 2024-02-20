# cuda and cublas utils

if(DCA_HAVE_GPU)
  set(srcs util_gpublas.cpp error_gpuBLAS.cpp info_gpu.cpp error_gpu.cpp)
  add_library(gpu_utils STATIC ${srcs})
  dca_gpu_runtime_link(gpu_utils)
  dca_gpu_blas_link(gpu_utils)
  target_link_libraries(gpu_utils PUBLIC magma::magma roc::hipblas) #"$<TARGET_PROPERTY:magma::magma,INTERFACE_LINK_LIBRARIES>")
  get_property(util_include_dirs TARGET gpu_utils PROPERTY INCLUDE_DIRECTORIES)
  message("gpu_utils includes: ${util_include_dirs}")
  get_property(interface_util_include_dirs TARGET gpu_utils PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
  message("gpu_utils interface includes: ${interface_util_include_dirs}")

endif()

