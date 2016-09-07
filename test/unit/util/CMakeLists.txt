# util unit tests

dca_add_gtest(ignore_unused_test
  GTEST_MAIN
  INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)

dca_add_gtest(print_type_test
  GTEST_MAIN
  INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)
