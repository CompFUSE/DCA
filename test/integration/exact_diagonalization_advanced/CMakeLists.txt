# Advanced exact diagonalization integration tests

dca_add_gtest(hamiltonian_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR};${PROJECT_SOURCE_DIR}
  LIBS function function_transform dca_io json time_and_frequency_domains cluster_domains enumerations quantum_domains ${LAPACK_LIBRARIES}
       ${DCA_GPU_LIBS} ${UTIL_LIBS} parallel_no_concurrency)

dca_add_gtest(ed_cluster_solver_single_site_test
  GTEST_MAIN
  CUDA
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS function function_transform json time_and_frequency_domains cluster_domains enumerations quantum_domains dca_io timer
       dca_algorithms ${LAPACK_LIBRARIES} ${DCA_THREADING_LIBS} ${DCA_GPU_LIBS} dca_hdf5 gpu_utils parallel_no_concurrency)

dca_add_gtest(ed_cluster_solver_four_site_test
  EXTENSIVE
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS function function_transform json time_and_frequency_domains cluster_domains enumerations quantum_domains dca_io timer
       dca_algorithms ${LAPACK_LIBRARIES} ${DCA_GPU_LIBS} ${DCA_THREADING_LIBS})
