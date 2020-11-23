dca_add_gtest(ringG_tp_accumulator_gpu_test
 FAST
 CUDA_MPI
 MPI MPI_NUMPROC ${DCA_TEST_GPU_COUNT}
 GTEST_MAIN
 INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
 LIBS     function ${DCA_LIBS} parallel_mpi_concurrency function)
