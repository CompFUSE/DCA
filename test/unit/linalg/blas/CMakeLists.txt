# linalg/util unit tests

if(DCA_HAVE_GPU)
  add_library(complex_op_test_kernels complex_op_test_kernels.cu)
  dca_gpu_device_link(complex_op_test_kernels)

  dca_add_gtest(complex_op_test
                GTEST_MAIN
                CUDA
		LIBS ${DCA_LIBS} complex_op_test_kernels)
endif()

