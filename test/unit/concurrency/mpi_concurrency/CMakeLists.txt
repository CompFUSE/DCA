# MPI concurrency's unit tests

dca_add_gtest(mpi_processor_grouping_test MPI MPI_NUMPROC 4)
dca_add_gtest(mpi_type_map_test GTEST_MAIN MPI MPI_NUMPROC 1)
