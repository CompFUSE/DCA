# Filename parameters' unit tests

dca_add_gtest(filename_parameters_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS json ${LAPACK_LIBRARIES})
