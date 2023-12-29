# util's unit tests

dca_add_gtest(ignore_test GTEST_MAIN)
dca_add_gtest(integer_division_test GTEST_MAIN)
dca_add_gtest(print_type_test GTEST_MAIN)
dca_add_gtest(plot_test
  INCLUDE_DIRS ${GNUPLOT_INTERFACE_INCLUDE_DIR}
  LIBS ${GNUPLOT_INTERFACE_LIBRARY})
dca_add_gtest(print_time_test GTEST_MAIN)
dca_add_gtest(timer_test GTEST_MAIN LIBS timer)
dca_add_gtest(pack_operations_test GTEST_MAIN)
dca_add_gtest(type_list_test GTEST_MAIN)
dca_add_gtest(type_utils_test GTEST_MAIN)
if(DCA_HAVE_GPU)
  dca_add_gtest(type_help_test GTEST_MAIN CUDA LIBS gpu_utils)
else()
  dca_add_gtest(type_help_test GTEST_MAIN)
endif()
add_subdirectory(containers)
