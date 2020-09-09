dca_add_gtest(dca_sp_DCA+_mpi_test
  EXTENSIVE
  MPI MPI_NUMPROC 8
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
  LIBS ${DCA_LIBS})

dca_add_gtest(dca_sp_DCA+_mpi_cuda_test
  EXTENSIVE
  MPI MPI_NUMPROC 8
  CUDA
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
  LIBS ${DCA_LIBS} ${DCA_CUDA_LIBS})

dca_add_gtest(dca_sp_DCA+_thread_test
  GTEST_MAIN
  EXTENSIVE
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
  LIBS ${DCA_LIBS})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/initial_sigma.hdf5
               ${CMAKE_CURRENT_BINARY_DIR}/initial_sigma.hdf5 COPYONLY)
