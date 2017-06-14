# Advanced exact diagonalization integration tests

dca_add_gtest(hamiltonian_test
  GTEST_MAIN
  LIBS function json time_and_frequency_domains cluster_domains quantum_domains ${LAPACK_LIBRARIES})
