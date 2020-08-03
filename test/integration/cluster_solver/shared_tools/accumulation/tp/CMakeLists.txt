dca_add_gtest(ringG_tp_accumulator_gpu_test
        EXTENSIVE
        GTEST_MAIN
        CUDA
        MPI MPI_NUMPROC 3
        INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
        LIBS     ${DCA_LIBS} parallel_mpi_concurrency
        )
