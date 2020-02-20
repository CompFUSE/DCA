# thread pool unit tests

dca_add_gtest(thread_pool_test GTEST_MAIN LIBS ${DCA_THREADING_LIBS})
dca_add_gtest(affinity_test GTEST_MAIN LIBS ${DCA_THREADING_LIBS})
