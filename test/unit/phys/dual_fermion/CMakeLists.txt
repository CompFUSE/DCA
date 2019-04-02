# Dual Fermion unit tests

dca_add_gtest(dual_self_energy_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function time_and_frequency_domains ${LAPACK_LIBRARIES})

dca_add_gtest(dual_to_lattice_self_energy_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function time_and_frequency_domains ${LAPACK_LIBRARIES})

dca_add_gtest(high_frequency_tails_test
  GTEST_MAIN
  LIBS function time_and_frequency_domains)

dca_add_gtest(high_frequency_tails_same_size_test
  GTEST_MAIN
  LIBS function time_and_frequency_domains)

dca_add_gtest(bare_dual_greens_function_test
  GTEST_MAIN
  LIBS function ${LAPACK_LIBRARIES})
