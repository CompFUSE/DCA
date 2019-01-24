# IO's unit tests

dca_add_gtest(buffer_test
  GTEST_MAIN)

dca_add_gtest(reader_test
  GTEST_MAIN
  LIBS function dca_hdf5 json ${HDF5_LIBRARIES})

dca_add_gtest(writer_test
  GTEST_MAIN
  LIBS function dca_hdf5 json ${HDF5_LIBRARIES})

dca_add_gtest(hdf5_reader_writer_test
  GTEST_MAIN
  LIBS dca_hdf5 ${HDF5_LIBRARIES})
