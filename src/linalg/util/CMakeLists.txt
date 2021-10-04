# cuda and cublas utils

if(DCA_HAVE_CUDA OR DCA_HAVE_HIP)
  set(srcs util_gpublas.cpp error_gpuBLAS.cpp info_gpu.cpp error_gpu.cpp)
  add_library(gpu_utils STATIC ${srcs})
#  target_include_directories(gpu_utils PUBLIC ${MAGMA_INCLUDE_DIR})
  target_link_libraries(gpu_utils PUBLIC magma::magma)
  if(DCA_HAVE_CUDA)
    target_link_libraries(gpu_utils PUBLIC CUDA::cudart CUDA::cublas)
    target_compile_definitions(gpu_utils PUBLIC DCA_HAVE_CUDA)
  endif()
  if(DCA_HAVE_HIP)
    target_link_libraries(gpu_utils PUBLIC hip::host)
    target_compile_definitions(gpu_utils PUBLIC DCA_HAVE_HIP)
  endif()
endif()

