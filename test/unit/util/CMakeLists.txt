# util's unit tests

dca_add_gtest(ignore_test GTEST_MAIN)

dca_add_gtest(print_type_test
  GTEST_MAIN
  INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)
