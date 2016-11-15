# Function library's unit tests

add_subdirectory(domains)

dca_add_gtest(function_library_test
  GTEST_MAIN
  LIBS function)
