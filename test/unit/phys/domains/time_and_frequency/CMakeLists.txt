# Time and frequency domain unit tests

dca_add_gtest(time_domain_test
  GTEST_MAIN
  LIBS time_and_frequency_domains)

dca_add_gtest(time_domain_left_oriented_test
  GTEST_MAIN
  LIBS time_and_frequency_domains)
