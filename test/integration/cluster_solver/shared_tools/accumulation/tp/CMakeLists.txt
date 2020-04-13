dca_add_gtest(ringG_tp_accumulator_gpu_test
        CUDA
        MPI MPI_NUMPROC 3
        GTEST_MAIN
        INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
        LIBS     ${DCA_LIBS} parallel_mpi_concurrency
        )
