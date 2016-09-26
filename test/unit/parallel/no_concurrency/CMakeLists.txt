# No-concurrency's unit tests

dca_add_gtest(no_concurrency_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS ${FFTW_LIBRARY})
dca_add_gtest(serial_collective_sum_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS ${FFTW_LIBRARY})
dca_add_gtest(serial_packing_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS ${FFTW_LIBRARY})
dca_add_gtest(serial_processor_grouping_test GTEST_MAIN)
