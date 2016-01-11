################################################################################
# test/DCA+/dca/mpi/CMakeLists.txt
################################################################################

# dca sp
file(COPY input.U=4_d=0.900_Nc=4A_T=1.0.DCA+_sp.json
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_gtest(NAME DCA+_sp_mpi_test
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
