# nfft integration tests

dca_add_gtest(dnfft_1d_test
  GTEST_MAIN
  FAST
  INCLUDE_DIRS ${FFTW_INCLUDE_DIRS}
  LIBS ${FFTW_LIBRARIES} time_and_frequency_domains random function nfft)
