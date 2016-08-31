# DCA(+) analysis system-level tests

if (DCA_WITH_MPI)
  # TODO: Avoid copying this file.
  file(COPY input_data.analysis_DCA+_mpi_test.hdf5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  add_gtest(NAME analysis_DCA+_mpi_test
  ADVANCED
  MPI MPI_NUMPROC 8
  INCLUDES ${DCA_INCLUDES}
  LIBS ${DCA_LIBRARIES})

target_compile_definitions(analysis_DCA+_mpi_test PRIVATE DCA_WITH_REDUCED_VERTEX_FUNCTION)
endif()
