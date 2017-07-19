# Function library's unit tests

add_subdirectory(domains)

dca_add_gtest(function_library_test
  GTEST_MAIN
  LIBS function)

dca_add_gtest(function_utils_test
    GTEST_MAIN
    LIBS function)

dca_add_gtest(set_to_zero_test GTEST_MAIN)
