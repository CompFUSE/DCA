# Dual Fermion unit tests

dca_add_gtest(high_frequency_tails_test
  GTEST_MAIN
  LIBS function time_and_frequency_domains)

dca_add_gtest(high_frequency_tails_same_size_test
  GTEST_MAIN
  LIBS function time_and_frequency_domains)

dca_add_gtest(compute_bare_dual_greens_function_test
  GTEST_MAIN
  LIBS function)
