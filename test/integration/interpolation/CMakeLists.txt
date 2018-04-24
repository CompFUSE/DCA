# Interpolation tests

dca_add_gtest(wannier_interpolation_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function ${LAPACK_LIBRARIES})
