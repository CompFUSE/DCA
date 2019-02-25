# LAPACK

if(DCA_HAVE_CUDA)
  CUDA_ADD_LIBRARY(lapack_kernels laset_gpu.cu multiply_diagonal_gpu.cu)
  target_compile_definitions(lapack_kernels PRIVATE DCA_HAVE_CUDA)
endif()

add_library(lapack silence_lapack.cpp)
target_include_directories(lapack PRIVATE ${DCA_ESSL_INCLUDE})
if(DCA_ESSL_INCLUDES)
  target_compile_definitions(lapack PUBLIC DCA_WITH_ESSL)
  target_include_directories(lapack PRIVATE ${DCA_ESSL_INCLUDES})
endif()
