# Function library's unit tests

add_subdirectory(domains)

dca_add_gtest(function_library_test
  GTEST_MAIN
  INCLUDE_DIRS ${PROJECT_SOURCE_DIR}
  LIBS function time_and_frequency_domains
  )
