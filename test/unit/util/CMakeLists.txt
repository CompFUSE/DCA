# util's unit tests

dca_add_gtest(ignore_test GTEST_MAIN)
dca_add_gtest(integer_division_test GTEST_MAIN)
dca_add_gtest(print_type_test GTEST_MAIN)
dca_add_gtest(plot_test
              INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
              LIBS gnuplot_interface)
