# Fourier transformation tests

dca_add_gtest(cluster_fourier_transform_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIRS}
  LIBS function ${LAPACK_LIBRARIES})
