# nfft integration tests

dca_add_gtest(dnfft_1d_test
  GTEST_MAIN
  FAST
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR};${PROJECT_SOURCE_DIR}
  LIBS FFTW::Double ${DCA_LIBS} time_and_frequency_domains random function nfft)
