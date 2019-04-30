dca_add_gtest(large_mpi_sum_test
              EXTENSIVE
              GTEST_MAIN
              LIBS parallel_mpi_concurrency function
              MPI MPI_NUMPROC 2)
