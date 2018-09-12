# util's unit tests

dca_add_gtest(call_once_per_loop_test GTEST_MAIN LIBS parallel_stdthread)
dca_add_gtest(ignore_test GTEST_MAIN)
dca_add_gtest(integer_division_test GTEST_MAIN)
dca_add_gtest(print_type_test GTEST_MAIN)
dca_add_gtest(plot_test
  INCLUDE_DIRS ${GNUPLOT_INTERFACE_INCLUDE_DIR}
  LIBS ${GNUPLOT_INTERFACE_LIBRARY})
dca_add_gtest(print_time_test GTEST_MAIN)
dca_add_gtest(timer_test GTEST_MAIN LIBS timer)
dca_add_gtest(pack_operations_test GTEST_MAIN)
