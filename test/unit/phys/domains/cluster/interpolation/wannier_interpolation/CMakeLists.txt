# Wannier interpolation unit tests
#
# Wannier interpolation is deprecated!

dca_add_gtest(wannier_interpolation_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NFFT_INCLUDE_DIR}
  LIBS json cluster_domains quantum_domains time_and_frequency_domains function ${NFFT_LIBRARY}
  ${FFTW_LIBRARY})
