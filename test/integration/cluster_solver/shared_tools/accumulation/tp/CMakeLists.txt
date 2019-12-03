dca_add_gtest(distributed_tp_accumulator_singleband_gpu_test
        CUDA
        GTEST_MAIN
        INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
        LIBS     ${DCA_LIBS}
        )

dca_add_gtest(distributed_tp_accumulator_gpu_test
  CUDA
  MPI MPI_NUMPROC 8
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS} parallel_mpi_concurrency
  )
