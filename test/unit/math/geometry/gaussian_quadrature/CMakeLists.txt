# gaussian quadrature tests

add_gtest(NAME compute_weights_and_abscissas_test
          GTEST_MAIN
          INCLUDES ${DCA_INCLUDES}
          LIBS gaussian_quadrature)
