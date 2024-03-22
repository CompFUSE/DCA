# No concurrency unit tests

set(optional_io_lib "")

if (DCA_HAVE_ADIOS2)
  set(optional_io_lib dca_adios2)
endif()

dca_add_gtest(no_concurrency_test
  GTEST_MAIN
  LIBS function
  ${optional_io_lib}
  ${DCA_LIBS}
  parallel_no_concurrency
  )

dca_add_gtest(serial_collective_sum_test
  GTEST_MAIN
  LIBS function)
dca_add_gtest(serial_packing_test
  GTEST_MAIN
  LIBS function)
dca_add_gtest(serial_processor_grouping_test GTEST_MAIN)
