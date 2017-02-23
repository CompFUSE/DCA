# math/util unit tests

dca_add_gtest(vector_operations_test
  GTEST_MAIN
  LIBS ${LAPACK_LIBRARIES})

dca_add_gtest(comparison_methods_test
  GTEST_MAIN
  LIBS ${LAPACK_LIBRARIES})
