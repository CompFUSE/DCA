# IO's unit tests

dca_add_gtest(reader_test
  GTEST_MAIN
  LIBS json_reader ${LAPACK_LIBRARIES})

dca_add_gtest(writer_test
  GTEST_MAIN
  LIBS json_writer ${LAPACK_LIBRARIES})
