# MPI concurrency unit tests

# If CUDA is enabled GPU >= the number of mpi proc must be available.

set(MPI_TEST_LIBS parallel_mpi_concurrency)
if(DCA_HAVE_GPU)
  set(MPI_TEST_LIBS "${MPI_TEST_LIBS};gpu_utils")
endif(DCA_HAVE_GPU)

dca_add_gtest(mpi_collective_max_test
  MPI MPI_NUMPROC 4
  LIBS ${MPI_TEST_LIBS})
dca_add_gtest(mpi_collective_min_test
  MPI MPI_NUMPROC 4
  LIBS ${MPI_TEST_LIBS})
dca_add_gtest(mpi_collective_sum_test
  MPI MPI_NUMPROC 6
  LIBS ${MPI_TEST_LIBS} function)
dca_add_gtest(mpi_packing_test
  MPI MPI_NUMPROC 1
  LIBS ${MPI_TEST_LIBS} function)
dca_add_gtest(mpi_processor_grouping_test
  MPI MPI_NUMPROC 4
  LIBS ${MPI_TEST_LIBS})
dca_add_gtest(mpi_type_map_test
  GTEST_MAIN
  MPI MPI_NUMPROC 1)
dca_add_gtest(mpi_concurrency_test
  MPI MPI_NUMPROC 4
  LIBS ${MPI_TEST_LIBS})
dca_add_gtest(mpi_gather_test
  MPI MPI_NUMPROC 6
  LIBS ${MPI_TEST_LIBS} function)
