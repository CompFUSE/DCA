################################################################################
# applications/dca/test/CMakeLists.txt
################################################################################

# dca sp DCA+ MPI
file(COPY input.dca_sp_DCA+_mpi_test.json
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY check_data.dca_sp_DCA+_mpi_test.hdf5
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

list(APPEND DCA_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}")

set(TEST_INCLUDES  "${PROJECT_SOURCE_DIR}/applications;${DCA_INCLUDES}") #use default types for test

add_gtest(NAME dca_sp_DCA+_mpi_test
  ADVANCED
  MPI MPI_NUMPROC 8
  INCLUDES "${TEST_INCLUDES}"
  LIBS     "${DCA_LIBRARIES}"
  )
