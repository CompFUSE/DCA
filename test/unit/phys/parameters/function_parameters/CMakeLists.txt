# Function parameters' unit tests

dca_add_gtest(function_parameters_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS ${LAPACK_LIBRARIES})
