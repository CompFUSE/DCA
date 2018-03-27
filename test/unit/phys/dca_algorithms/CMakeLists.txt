# DCA algorithms unit tests

dca_add_gtest(compute_band_structure_test
  GTEST_MAIN
  LIBS dca_algorithms json function cluster_domains time_and_frequency_domains enumerations
       quantum_domains
  ${LAPACK_LIBRARIES})

dca_add_gtest(compute_free_greens_function_test
  GTEST_MAIN
  LIBS function ${LAPACK_LIBRARIES})

dca_add_gtest(compute_greens_function_test
  GTEST_MAIN
  LIBS function ${LAPACK_LIBRARIES})
