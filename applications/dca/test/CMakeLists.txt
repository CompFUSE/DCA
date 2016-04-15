################################################################################
# applications/dca/test/CMakeLists.txt
################################################################################
set(TEST_INCLUDES
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${DCA_INCLUDES}
        )


set(TEST_INCLUDES  "${CMAKE_CURRENT_SOURCE_DIR};${DCA_INCLUDES}") 

add_gtest(NAME dca_sp_DCA+_mpi_test
  ADVANCED
  MPI MPI_NUMPROC 8
  INCLUDES "${TEST_INCLUDES}"
  LIBS     "${DCA_LIBRARIES}"
  )

add_gtest(NAME dca_sp_DCA+_pthread_test
  ADVANCED
  MPI MPI_NUMPROC 1
  INCLUDES "${TEST_INCLUDES}"
  LIBS     "${DCA_LIBRARIES}"
  )
