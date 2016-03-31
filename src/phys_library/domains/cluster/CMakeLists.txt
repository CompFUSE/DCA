set(INCLUDES ${PROJECT_SOURCE_DIR}/src)
# message("INCLUDE = ${INCLUDES}")
set(LIBS ${LAPACK_LIBRARIES})
# message("LIBS = ${LIBS}")

add_gtest(
  NAME cluster_operations_test
  GTEST_MAIN
  INCLUDES "${INCLUDES}"
  LIBS "${LIBS}"
  )
