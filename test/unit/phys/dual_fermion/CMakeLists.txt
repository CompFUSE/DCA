# Dual Fermion unit tests

dca_add_gtest(dual_self_energy_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function time_and_frequency_domains ${LAPACK_LIBRARIES})

dca_add_gtest(dual_to_lattice_self_energy_test
  GTEST_MAIN
  LIBS function time_and_frequency_domains ${LAPACK_LIBRARIES})
