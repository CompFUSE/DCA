# domains/cluster unit tests

dca_add_gtest(cluster_operations_test
  GTEST_MAIN
  LIBS ${LAPACK_LIBRARIES} ${DCA_CUDA_LIBS})

dca_add_gtest(momentum_exchage_domain_test
  GTEST_MAIN
  LIBS cluster_domains)

# deprecated (requires NFFT)
# add_subdirectory(interpolation/wannier_interpolation)
