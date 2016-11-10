# Double-counting parameters' unit tests

dca_add_gtest(double_counting_parameters_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS json ${LAPACK_LIBRARIES})
