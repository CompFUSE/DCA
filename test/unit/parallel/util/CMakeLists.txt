# parallel/util unit tests

dca_add_gtest(call_once_per_loop_test GTEST_MAIN LIBS ${DCA_THREADING_LIBS})

dca_add_gtest(get_bounds_test
  GTEST_MAIN
  LIBS function parallel_util)
