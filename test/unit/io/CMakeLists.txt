# IO's unit tests

dca_add_gtest(buffer_test
  GTEST_MAIN)

dca_add_gtest(reader_test
  MPI
  LIBS function ${DCA_LIBS} dca_io)

dca_add_gtest(writer_test
  MPI
  LIBS function dca_io parallel_mpi_concurrency parallel_no_concurrency)

dca_add_gtest(json_reader_test
  GTEST_MAIN
  LIBS json)

dca_add_gtest(reader_writer_test
  LIBS dca_io parallel_mpi_concurrency function ${DCA_THREADING_LIBS})

dca_add_gtest(io_types_test
  GTEST_MAIN
  LIBS dca_io)

if (DCA_HAVE_ADIOS2)
  if (DCA_HAVE_MPI)
    dca_add_gtest(adios2_reader_writer_test
      LIBS parallel_mpi_concurrency dca_io function)
    dca_add_gtest(adios2_parallelIO_test
      MPI MPI_NUMPROC 6
      LIBS dca_io parallel_mpi_concurrency function parallel_util)
  else()
    dca_add_gtest(adios2_reader_writer_test
      LIBS parallel_no_concurrency dca_io function)
  endif()
endif()

