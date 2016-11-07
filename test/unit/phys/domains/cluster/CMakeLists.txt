dca_add_gtest(cluster_operations_test
  GTEST_MAIN
  LIBS ${LAPACK_LIBRARIES})

add_subdirectory(interpolation/wannier_interpolation)
