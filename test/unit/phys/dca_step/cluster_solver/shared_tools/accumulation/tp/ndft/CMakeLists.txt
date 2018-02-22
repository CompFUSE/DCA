# test/unit/phys/dca_step/cluster_solver/shared_tools/tp/ndft

dca_add_gtest(cached_ndft_test
  GTEST_MAIN
  FAST
  LIBS ${LAPACK_LIBRARIES} time_and_frequency_domains random function)

dca_add_gtest(cached_ndft_gpu_test
  GTEST_MAIN
  CUDA
  FAST
  LIBS ${DCA_CUDA_LIBS} time_and_frequency_domains random function mc_tools_kernels)
