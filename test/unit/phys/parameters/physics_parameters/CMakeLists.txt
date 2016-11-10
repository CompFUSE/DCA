# Physics parameters' unit tests

dca_add_gtest(physics_parameters_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS json ${LAPACK_LIBRARIES})
