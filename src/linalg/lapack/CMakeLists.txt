# LAPACK

if(DCA_HAVE_CUDA)
  CUDA_ADD_LIBRARY(lapack_kernels laset_gpu.cu multiply_diagonal_gpu.cu)
  target_compile_definitions(lapack_kernels PRIVATE -DDCA_HAVE_CUDA)
endif()
