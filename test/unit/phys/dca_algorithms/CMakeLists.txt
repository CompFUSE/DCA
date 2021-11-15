# DCA algorithms unit tests

dca_add_gtest(compute_band_structure_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS dca_algorithms json function cluster_domains enumerations time_and_frequency_domains quantum_domains ${UTIL_LIBS}
  ${LAPACK_LIBRARIES} ${DCA_GPU_LIBS})

# \todo I think this should be able to run without GPU
dca_add_gtest(compute_free_greens_function_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  CUDA
  LIBS function ${DCA_THREADING_LIBS} ${LAPACK_LIBRARIES} ${UTIL_LIBS} ${DCA_GPU_LIBS} parallel_util)

dca_add_gtest(compute_greens_function_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function ${DCA_THREADING_LIBS} ${LAPACK_LIBRARIES} ${UTIL_LIBS} ${DCA_GPU_LIBS} parallel_util)
