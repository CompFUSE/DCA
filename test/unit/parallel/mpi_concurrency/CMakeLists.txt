# MPI concurrency unit tests

dca_add_gtest(mpi_collective_max_test
  MPI MPI_NUMPROC 4)
dca_add_gtest(mpi_collective_min_test
  MPI MPI_NUMPROC 4)
dca_add_gtest(mpi_collective_sum_test
  MPI MPI_NUMPROC 8
  LIBS function ${DCA_CUDA_LIBS})
dca_add_gtest(mpi_collective_gather_test
  MPI MPI_NUMPROC 4
  LIBS function parallel_util)
dca_add_gtest(mpi_packing_test
  MPI MPI_NUMPROC 1
  LIBS function)
dca_add_gtest(mpi_processor_grouping_test
  MPI MPI_NUMPROC 4)
dca_add_gtest(mpi_type_map_test
  GTEST_MAIN
  MPI MPI_NUMPROC 1)
dca_add_gtest(mpi_concurrency_test
  MPI MPI_NUMPROC 4)
