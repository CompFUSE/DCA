# IO's unit tests

dca_add_gtest(buffer_test
  GTEST_MAIN)

dca_add_gtest(reader_test
  GTEST_MAIN
  LIBS function dca_hdf5 json ${HDF5_LIBRARIES})

dca_add_gtest(writer_test
  GTEST_MAIN
  LIBS function dca_hdf5 json ${HDF5_LIBRARIES})

dca_add_gtest(json_reader_test
        GTEST_MAIN
        LIBS json)

dca_add_gtest(reader_writer_test
  GTEST_MAIN
  LIBS dca_hdf5;json;function)
