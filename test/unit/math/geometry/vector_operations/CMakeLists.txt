dca_add_gtest(elementary_vector_operations_test
  GTEST_MAIN
  INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/src
  LIBS ${LAPACK_LIBRARIES})
