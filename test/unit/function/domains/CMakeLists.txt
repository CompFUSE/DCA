# Unit tests for domain classes

dca_add_gtest(dmn_test
  GTEST_MAIN
  INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)

dca_add_gtest(reduced_domain_test
  GTEST_MAIN
  INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)

dca_add_gtest(dmn_0_test
  GTEST_MAIN
  INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)
