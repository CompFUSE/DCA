# cuda and cublas utils

if(DCA_HAVE_CUDA)
  set(srcs util_cublas.cpp error_cublas.cpp info_cuda.cpp error_cuda.cpp)
  add_library(cuda_utils STATIC ${srcs})
  target_include_directories(cuda_utils PUBLIC ${MAGMA_INCLUDE_DIR};${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
  target_compile_definitions(cuda_utils PUBLIC DCA_HAVE_CUDA)
else()
  add_library(cuda_utils STATIC util_cublas.cpp)
endif()

