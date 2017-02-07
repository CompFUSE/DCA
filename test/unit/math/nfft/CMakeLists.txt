# math/nfft unit tests

dca_add_gtest(dnfft_1d_test
  GTEST_MAIN
  LIBS ${FFTW_LIBRARY};time_and_frequency_domains;random;function
  )
