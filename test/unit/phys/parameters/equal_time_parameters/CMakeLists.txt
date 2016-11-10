# Equal-time parameters' unit tests

dca_add_gtest(equal_time_parameters_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS json ${LAPACK_LIBRARIES})
