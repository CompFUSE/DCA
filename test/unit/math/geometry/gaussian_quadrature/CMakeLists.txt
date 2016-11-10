# Gaussian quadrature unit tests

dca_add_gtest(compute_weights_and_abscissas_test
  GTEST_MAIN
  LIBS gaussian_quadrature tetrahedron_mesh ${LAPACK_LIBRARIES} ${SIMPLEX_GM_RULE_LIBRARY})
