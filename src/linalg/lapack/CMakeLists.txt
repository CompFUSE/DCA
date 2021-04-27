# LAPACK

if(DCA_HAVE_CUDA)
  add_library(lapack_kernels laset_gpu.cu multiply_diagonal_gpu.cu)
  set_target_properties( lapack_kernels
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

                     target_compile_definitions(lapack_kernels PRIVATE DCA_HAVE_CUDA)
                     target_link_libraries(lapack_kernels PRIVATE CUDA::cublas)
endif()

add_library(lapack silence_lapack.cpp)
target_include_directories(lapack PRIVATE ${DCA_ESSL_INCLUDE})
if(DCA_ESSL_INCLUDES)
  target_compile_definitions(lapack PUBLIC DCA_WITH_ESSL)
  target_include_directories(lapack PRIVATE ${DCA_ESSL_INCLUDES})
endif()
