# gaussian quadrature unit tests

add_gtest(NAME compute_weights_and_abscissas_test
          GTEST_MAIN
          LIBS "gaussian_quadrature;${LAPACK_LIBRARIES}")
