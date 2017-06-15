# Advanced exact diagonalization integration tests

dca_add_gtest(hamiltonian_test
  GTEST_MAIN
  LIBS function json time_and_frequency_domains cluster_domains quantum_domains ${LAPACK_LIBRARIES})

dca_add_gtest(ed_cluster_solver_single_site_test
  GTEST_MAIN
  INCLUDES
  LIBS function json time_and_frequency_domains cluster_domains quantum_domains hdf5 timer
       dca_algorithms ${LAPACK_LIBRARIES})

dca_add_gtest(ed_cluster_solver_four_site_test
  EXTENSIVE
  GTEST_MAIN
  INCLUDES
  LIBS function json time_and_frequency_domains cluster_domains quantum_domains hdf5 timer
       dca_algorithms ${LAPACK_LIBRARIES})
