# test/unit/phys/dca_step/cluster_solver/shared_tools/tp/ndft

dca_add_gtest(cached_ndft_cpu_test
  GTEST_MAIN
  FAST
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS ${LAPACK_LIBRARIES} time_and_frequency_domains random function)

dca_add_gtest(cached_ndft_gpu_test
  GTEST_MAIN
  CUDA
  FAST
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS ${DCA_CUDA_LIBS} time_and_frequency_domains random function mc_kernels cuda_utils)
# ${LAPACK_LIBRARIES}
