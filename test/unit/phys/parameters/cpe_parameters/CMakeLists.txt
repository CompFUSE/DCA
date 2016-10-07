# Continuous-pole-expansion (CPE) parameters' unit tests

dca_add_gtest(cpe_parameters_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS json_reader ${LAPACK_LIBRARIES})
