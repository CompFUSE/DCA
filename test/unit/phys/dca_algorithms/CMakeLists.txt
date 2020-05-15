# DCA algorithms unit tests

dca_add_gtest(compute_band_structure_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS dca_algorithms json function cluster_domains enumerations time_and_frequency_domains quantum_domains
  ${LAPACK_LIBRARIES} ${DCA_CUDA_LIBS} cuda_utils)

dca_add_gtest(compute_free_greens_function_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function ${DCA_THREADING_LIBS} parallel_util ${LAPACK_LIBRARIES} ${DCA_CUDA_LIBS} cuda_utils)

dca_add_gtest(compute_greens_function_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function ${DCA_THREADING_LIBS} parallel_util ${LAPACK_LIBRARIES} ${DCA_CUDA_LIBS} cuda_utils)
