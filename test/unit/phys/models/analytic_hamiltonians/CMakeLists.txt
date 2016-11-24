# Analytic Hamiltonians unit tests

dca_add_gtest(bilayer_lattice_test
  GTEST_MAIN
  LIBS function)

dca_add_gtest(square_lattice_test
  GTEST_MAIN
  LIBS function ${LAPACK_LIBRARIES})

dca_add_gtest(triangular_lattice_test
  GTEST_MAIN
  LIBS function)
