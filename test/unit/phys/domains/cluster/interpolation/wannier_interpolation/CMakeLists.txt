# Wannier interpolation unit tests

dca_add_gtest(wannier_interpolation_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NFFT_INCLUDE_DIR}
  LIBS json_reader cluster_domains quantum_domains time_and_frequency_domains function
  ${NFFT_LIBRARY} ${FFTW_LIBRARY})
