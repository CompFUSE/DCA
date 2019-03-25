# DCA + DF integration tests

dca_add_gtest(non_translational_invariant_dispersion_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function ${LAPACK_LIBRARIES})
