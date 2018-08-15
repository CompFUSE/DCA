# cuda and cublas utils

if(DCA_HAVE_CUDA)
  add_library(cuda_utils STATIC util_cublas.cpp error_cublas.cpp info_cuda.cpp error_cuda.cpp)
  target_include_directories(cuda_utils PUBLIC ${MAGMA_INCLUDE_DIR};${CUDA_TOOLKIT_INCLUDE})
  target_compile_definitions(cuda_utils PUBLIC DCA_HAVE_CUDA)
endif()
