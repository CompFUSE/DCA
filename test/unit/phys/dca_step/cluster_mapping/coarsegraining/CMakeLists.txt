# Coarsegraining unit tests

dca_add_gtest(coarsegraining_sp_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS json function cluster_domains time_and_frequency_domains quantum_domains gaussian_quadrature
       tetrahedron_mesh coarsegraining enumerations parallel_util parallel_stdthread lapack
       ${LAPACK_LIBRARIES})
