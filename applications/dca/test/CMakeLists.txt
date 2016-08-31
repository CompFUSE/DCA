# DCA(+) calculation system-level tests

if (DCA_HAVE_MPI)
  add_gtest(NAME dca_sp_DCA+_mpi_test
    ADVANCED
    MPI MPI_NUMPROC 8
    INCLUDES ${DCA_INCLUDES}
    LIBS ${DCA_LIBRARIES})
  target_compile_definitions(dca_sp_DCA+_mpi_test PRIVATE DCA_WITH_REDUCED_VERTEX_FUNCTION)

  if (DCA_HAVE_CUDA)
    add_gtest(NAME dca_sp_DCA+_mpi_cuda_test
      ADVANCED
      MPI MPI_NUMPROC 8
      CUDA
      INCLUDES ${DCA_INCLUDES}
      LIBS ${DCA_LIBRARIES})
    target_compile_definitions(dca_sp_DCA+_mpi_cuda_test PRIVATE DCA_WITH_REDUCED_VERTEX_FUNCTION)
  endif()
endif()

if (DCA_HAVE_PTHREADS)
  add_gtest(NAME dca_sp_DCA+_pthread_test
    ADVANCED
    MPI MPI_NUMPROC 1
    INCLUDES ${DCA_INCLUDES}
    LIBS ${DCA_LIBRARIES})
  target_compile_definitions(dca_sp_DCA+_pthread_test PRIVATE DCA_WITH_REDUCED_VERTEX_FUNCTION)
endif()
