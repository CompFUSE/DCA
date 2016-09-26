# DCA(+) analysis system-level tests

# TODO: Avoid copying this file.
file(COPY input_data.analysis_DCA+_mpi_test.hdf5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

dca_add_gtest(analysis_DCA+_mpi_test
  EXTENSIVE
  MPI MPI_NUMPROC 8
  PTHREADS
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
  LIBS ${DCA_LIBS})
