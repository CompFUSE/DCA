################################################################################
# applications/analysis/test/CMakeLists.txt
################################################################################

file(COPY input.analysis_DCA+_mpi_test.json
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY input_data.analysis_DCA+_mpi_test.hdf5
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY check_data.analysis_DCA+_mpi_test.hdf5
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


set(TEST_INCLUDES  "${CMAKE_CURRENT_SOURCE_DIR};${DCA_INCLUDES}")

add_gtest(NAME analysis_DCA+_mpi_test
  ADVANCED
  MPI MPI_NUMPROC 8
  INCLUDES "${TEST_INCLUDES}"
  LIBS     "${DCA_LIBRARIES}"
  )
