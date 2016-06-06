set(INCLUDES ${PROJECT_SOURCE_DIR}/src)
set(LIBS ${LAPACK_LIBRARIES})

add_gtest(
  NAME cluster_operations_test
  GTEST_MAIN
  INCLUDES "${INCLUDES}"
  LIBS "${LIBS}")
