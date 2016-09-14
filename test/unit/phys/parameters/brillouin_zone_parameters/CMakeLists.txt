# Brillouin zone parameters' unit tests

dca_add_gtest(brillouin_zone_parameters_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS ${LAPACK_LIBRARIES})
