# linalg/util unit tests

dca_add_gtest(check_RC_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS} ${DCA_GPU_LIBS})

dca_add_gtest(magma_batched_gemm_test
	      GTEST_MAIN
              CUDA
              LIBS lapack gpu_utils magma::magma magma::sparse ${DCA_EXTERNAL_LIBS})

dca_add_gtest(stream_handle_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS})

dca_add_gtest(pinned_allocator_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS} ${DCA_GPU_LIBS})

dca_add_gtest(copy_test
  GTEST_MAIN
  CUDA
  LIBS ${DCA_LIBS})

if(DCA_HAVE_CUDA)
dca_add_gtest(complex_op_cuda_test
  GTEST_MAIN
  CUDA
  LIBS ${DCA_GPU_LIBS})

dca_gpu_runtime_link(complex_op_cuda_test)
endif()
