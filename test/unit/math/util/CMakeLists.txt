# math/util unit tests
dca_add_gtest(phase_test
  GTEST_MAIN)

dca_add_gtest(vector_operations_test
  GTEST_MAIN
  LIBS ${MAGMA_LIBRARY} ${LAPACK_LIBRARIES} ${DCA_GPU_LIBS})

dca_add_gtest(comparison_methods_test
  GTEST_MAIN
  LIBS ${LAPACK_LIBRARIES})
