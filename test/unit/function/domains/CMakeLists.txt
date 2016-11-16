# Unit tests for the domain classes

dca_add_gtest(dmn_test GTEST_MAIN)
dca_add_gtest(dmn_0_test
  GTEST_MAIN
  LIBS function)
dca_add_gtest(domain_test
  GTEST_MAIN
  LIBS function)
dca_add_gtest(reduced_domain_test GTEST_MAIN)
dca_add_gtest(dmn_variadic_test
  GTEST_MAIN
  LIBS function)
