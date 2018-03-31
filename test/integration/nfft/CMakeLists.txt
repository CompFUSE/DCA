# nfft integration tests

dca_add_gtest(dnfft_1d_test
  GTEST_MAIN
  FAST
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS ${FFTW_LIBRARY} time_and_frequency_domains random function nfft)
