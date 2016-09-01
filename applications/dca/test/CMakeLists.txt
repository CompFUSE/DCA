# DCA(+) calculation system-level tests

dca_add_gtest(dca_sp_DCA+_mpi_test
  EXTENSIVE
  MPI MPI_NUMPROC 8
  INCLUDE_DIRS ${DCA_INCLUDES}
  LIBS ${DCA_LIBRARIES})

dca_add_gtest(dca_sp_DCA+_mpi_cuda_test
  EXTENSIVE
  MPI MPI_NUMPROC 8
  CUDA
  INCLUDE_DIRS ${DCA_INCLUDES}
  LIBS ${DCA_LIBRARIES})

dca_add_gtest(dca_sp_DCA+_pthread_test
  EXTENSIVE
  MPI MPI_NUMPROC 1
  PTHREADS
  INCLUDE_DIRS ${DCA_INCLUDES}
  LIBS ${DCA_LIBRARIES})
