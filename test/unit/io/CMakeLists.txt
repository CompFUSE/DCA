# IO's unit tests

dca_add_gtest(reader_test
  GTEST_MAIN
  LIBS hdf5 json ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})

dca_add_gtest(writer_test
  GTEST_MAIN
  LIBS hdf5 json ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})
