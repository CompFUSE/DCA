################################################################################
# applications/dca/test/CMakeLists.txt
################################################################################

# dca sp DCA+ MPI
file(COPY input.dca_sp_DCA+_mpi_test.json
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY check_data.dca_sp_DCA+_mpi_test.hdf5
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

list(APPEND DCA_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}")

add_gtest(NAME dca_sp_DCA+_mpi_test
  ADVANCED
  MPI MPI_NUMPROC 8
  INCLUDES "${DCA_INCLUDES}"
  LIBS     "${DCA_LIBRARIES}"
  )

# dca tp
# file(COPY input.U=4_d=0.900_Nc=4A_T=1.0.DCA+_tp.json
#   DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# file(COPY test_data.U=4_d=0.900_Nc=4A_T=1.0.DCA+_sp.hdf5
#   DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# add_gtest(NAME DCA+_tp_mpi_test
#   ADVANCED
#   MPI MPI_NUMPROC 8
#   INCLUDES "${DCA_INCLUDES}"
#   LIBS     "${DCA_LIBRARIES}"
#   )
