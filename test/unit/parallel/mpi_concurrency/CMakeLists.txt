# MPI concurrency unit tests

if (DCA_HAVE_CUDA)
  set(TEST_LIBS "parallel_mpi_concurrency;cuda_utils;${CUDA_CUDART_LIBRARY}")
else()
  set(TEST_LIBS parallel_mpi_concurrency)
endif()

dca_add_gtest(mpi_collective_max_test
  MPI MPI_NUMPROC 4
  LIBS ${TEST_LIBS})
dca_add_gtest(mpi_collective_min_test
  MPI MPI_NUMPROC 4
  LIBS ${TEST_LIBS})
dca_add_gtest(mpi_collective_sum_test
  MPI MPI_NUMPROC 8
  LIBS ${TEST_LIBS};function)
dca_add_gtest(mpi_packing_test
  MPI MPI_NUMPROC 1
  LIBS ${TEST_LIBS} function)
dca_add_gtest(mpi_processor_grouping_test
  MPI MPI_NUMPROC 4
  LIBS ${TEST_LIBS})
dca_add_gtest(mpi_type_map_test
  GTEST_MAIN
  MPI MPI_NUMPROC 1)
dca_add_gtest(mpi_concurrency_test
  MPI MPI_NUMPROC 4
  LIBS ${TEST_LIBS})
