# cuda and cublas utils

if(DCA_HAVE_CUDA)
  add_library(cuda_utils STATIC util_cublas.cpp error_cublas.cpp info_cuda.cpp error_cuda.cpp)
endif()
