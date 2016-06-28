################################################################################
# applications/dca/test/CMakeLists.txt
################################################################################

if (DCA_WITH_MPI)
  add_gtest(NAME dca_sp_DCA+_mpi_test
            ADVANCED
            MPI MPI_NUMPROC 8
            INCLUDES ${DCA_INCLUDES}
            LIBS ${DCA_LIBRARIES})

  if (DCA_WITH_CUDA)
    add_gtest(NAME dca_sp_DCA+_mpi_cuda_test
              ADVANCED
              MPI MPI_NUMPROC 8
              CUDA
              INCLUDES ${DCA_INCLUDES}
              LIBS ${DCA_LIBRARIES})
  endif()
endif()

if (DCA_WITH_PTHREADS)
  add_gtest(NAME dca_sp_DCA+_pthread_test
            ADVANCED
            MPI MPI_NUMPROC 1
            INCLUDES ${DCA_INCLUDES}
            LIBS ${DCA_LIBRARIES})
endif()
