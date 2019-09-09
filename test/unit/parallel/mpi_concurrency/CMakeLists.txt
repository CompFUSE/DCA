# MPI concurrency unit tests

dca_add_gtest(mpi_collective_max_test
  MPI MPI_NUMPROC 4
  LIBS parallel_mpi_concurrency)
dca_add_gtest(mpi_collective_min_test
  MPI MPI_NUMPROC 4
  LIBS parallel_mpi_concurrency)
dca_add_gtest(mpi_collective_sum_test
  MPI MPI_NUMPROC 8
  LIBS parallel_mpi_concurrency;function)
dca_add_gtest(mpi_packing_test
  MPI MPI_NUMPROC 1
  LIBS parallel_mpi_concurrency function)
dca_add_gtest(mpi_processor_grouping_test
  MPI MPI_NUMPROC 4
  LIBS parallel_mpi_concurrency)
dca_add_gtest(mpi_type_map_test
  GTEST_MAIN
  MPI MPI_NUMPROC 1)
dca_add_gtest(mpi_concurrency_test
  MPI MPI_NUMPROC 4
  LIBS parallel_mpi_concurrency)
dca_add_gtest(mpi_gather_test
  MPI MPI_NUMPROC 8
  LIBS parallel_mpi_concurrency function)
