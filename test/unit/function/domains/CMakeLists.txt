# Unit tests for domain classes

add_gtest(NAME dmn_test
  GTEST_MAIN
  INCLUDES "${PROJECT_SOURCE_DIR}/include")

add_gtest(NAME reduced_domain_test
  GTEST_MAIN
  INCLUDES "${PROJECT_SOURCE_DIR}/include")
