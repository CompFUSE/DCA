# MPI-concurrency's unit tests

dca_add_gtest(mpi_collective_max_test MPI MPI_NUMPROC 4)
dca_add_gtest(mpi_collective_min_test MPI MPI_NUMPROC 4)
dca_add_gtest(mpi_collective_sum_test
  MPI MPI_NUMPROC 8
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function ${FFTW_LIBRARY})
dca_add_gtest(mpi_packing_test
  MPI MPI_NUMPROC 4
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS ${FFTW_LIBRARY})
dca_add_gtest(mpi_processor_grouping_test MPI MPI_NUMPROC 4)
dca_add_gtest(mpi_type_map_test GTEST_MAIN MPI MPI_NUMPROC 1)
dca_add_gtest(mpi_concurrency_test
  MPI MPI_NUMPROC 4
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS ${FFTW_LIBRARY})
