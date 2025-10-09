# Fourier transformation tests

dca_add_gtest(cluster_fourier_transform_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function function_transform ${LAPACK_LIBRARIES} ${UTIL_LIBS})
