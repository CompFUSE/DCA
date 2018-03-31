# (Phys) deconvolution unit tests

dca_add_gtest(deconvolution_routines_test
  GTEST_MAIN
  INCLUDE_DIRS ${SIMPLEX_GM_RULE_INCLUDE_DIR}
  LIBS json function cluster_domains time_and_frequency_domains quantum_domains gaussian_quadrature
       tetrahedron_mesh coarsegraining ${LAPACK_LIBRARIES})
