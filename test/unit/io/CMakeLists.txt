# IO's unit tests

dca_add_gtest(reader_test
  GTEST_MAIN
  LIBS function dca_hdf5 json ${HDF5_LIBRARIES})

dca_add_gtest(writer_test
  GTEST_MAIN
  MPI MPI_NUMPROC 1
  LIBS function dca_hdf5 json ${HDF5_LIBRARIES})
