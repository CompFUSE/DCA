# Dual Fermion unit tests

dca_add_gtest(dual_self_energy_test
  GTEST_MAIN
  LIBS function time_and_frequency_domains ${LAPACK_LIBRARIES})
