add_subdirectory(mpi_concurrency)

dca_add_gtest(collective_sum_interface_mpi_test
  MPI MPI_NUMPROC 8
  INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src ${FFTW_INCLUDE_DIR}
  LIBS ${FFTW_LIBRARY})
