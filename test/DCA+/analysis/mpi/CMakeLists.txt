################################################################################
# test/DCA+/analysis/mpi/CMakeLists.txt
################################################################################

file(COPY input.U=4_d=0.900_Nc=4A_T=1.0.DCA+_analysis_ppSC.json
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY test_data.U=4_d=0.900_Nc=4A_T=1.0.DCA+_tp.hdf5
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_gtest(NAME DCA+_analysis_ppSC_mpi_test
  ADVANCED
  MPI MPI_NUMPROC 8
  INCLUDES "${DCA_INCLUDES}"
  LIBS     "${DCA_LIBRARIES}"
  )
