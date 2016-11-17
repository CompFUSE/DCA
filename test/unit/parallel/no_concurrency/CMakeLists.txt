# No concurrency unit tests

dca_add_gtest(no_concurrency_test
  GTEST_MAIN
  LIBS function)
dca_add_gtest(serial_collective_sum_test
  GTEST_MAIN
  LIBS function)
dca_add_gtest(serial_packing_test
  GTEST_MAIN
  LIBS function)
dca_add_gtest(serial_processor_grouping_test GTEST_MAIN)
