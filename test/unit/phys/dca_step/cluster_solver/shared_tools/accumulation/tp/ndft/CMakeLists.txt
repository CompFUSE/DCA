# test/unit/phys/dca_step/cluster_solver/shared_tools/tp/ndft

dca_add_gtest(cached_ndft_test
  GTEST_MAIN
  FAST
  LIBS ${DCA_LIBS})

dca_add_gtest(cached_ndft_gpu_test
  GTEST_MAIN
  CUDA
  FAST
  LIBS ${DCA_LIBS})
