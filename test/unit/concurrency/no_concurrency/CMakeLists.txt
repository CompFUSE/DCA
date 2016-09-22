# No-concurrency's unit tests

dca_add_gtest(serial_processor_grouping_test GTEST_MAIN)
dca_add_gtest(serial_collective_sum_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS ${FFTW_LIBRARY})
