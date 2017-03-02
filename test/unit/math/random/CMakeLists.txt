# Random number library unit tests

dca_add_gtest(random_test
  GTEST_MAIN
  LIBS random)

dca_add_gtest(rng_vector_construction_test
  GTEST_MAIN
  LIBS random)

# Random number library utility functions
dca_add_gtest(random_utils_test
  GTEST_MAIN
  LIBS random)

dca_add_gtest(random_utils_exceptions_test
  GTEST_MAIN
  LIBS random)

# Standard random number library
dca_add_gtest(std_random_wrapper_default_seed_test
  GTEST_MAIN
  LIBS random)

dca_add_gtest(std_random_wrapper_custom_seed_test
  GTEST_MAIN
  LIBS random)

dca_add_gtest(std_random_wrapper_unique_seeds_test
  GTEST_MAIN
  LIBS random)
