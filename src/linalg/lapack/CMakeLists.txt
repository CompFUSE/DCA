# LAPACK

if(DCA_HAVE_GPU)
  add_library(lapack_kernels laset_gpu.cu multiply_diagonal_gpu.cu)
  dca_gpu_device_link(lapack_kernels)
  dca_gpu_blas_link(lapack_kernels)
  dca_gpu_runtime_link(lapack_kernels)
endif()

add_library(lapack silence_lapack.cpp)
target_include_directories(lapack PRIVATE ${DCA_ESSL_INCLUDE})
if(DCA_ESSL_INCLUDES)
  target_compile_definitions(lapack PUBLIC DCA_WITH_ESSL)
  target_include_directories(lapack PRIVATE ${DCA_ESSL_INCLUDES})
endif()
