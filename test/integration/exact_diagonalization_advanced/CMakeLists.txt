# Advanced exact diagonalization integration tests

dca_add_gtest(hamiltonian_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIRS}
  LIBS function json time_and_frequency_domains cluster_domains enumerations quantum_domains ${LAPACK_LIBRARIES}
       ${DCA_CUDA_LIBS} cuda_utils)

dca_add_gtest(ed_cluster_solver_single_site_test
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
  LIBS function json time_and_frequency_domains cluster_domains enumerations quantum_domains dca_hdf5 timer
  dca_algorithms ${LAPACK_LIBRARIES} ${DCA_CUDA_LIBS} ${DCA_THREADING_LIBS} cuda_utils)

dca_add_gtest(ed_cluster_solver_four_site_test
  EXTENSIVE
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS} 
  LIBS function json time_and_frequency_domains cluster_domains enumerations quantum_domains dca_hdf5 timer
       dca_algorithms ${LAPACK_LIBRARIES} ${DCA_CUDA_LIBS} ${DCA_THREADING_LIBS} cuda_utils)
