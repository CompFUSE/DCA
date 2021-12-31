# IO's unit tests

dca_add_gtest(buffer_test
  GTEST_MAIN)

dca_add_gtest(reader_test
  GTEST_MAIN
  LIBS function dca_io)

dca_add_gtest(writer_test
  GTEST_MAIN
  LIBS function dca_io)

dca_add_gtest(json_reader_test
        GTEST_MAIN
        LIBS json)

dca_add_gtest(reader_writer_test
  LIBS dca_io function ${DCA_THREADING_LIBS})

if (DCA_HAVE_ADIOS2)
  dca_add_gtest(adios2_reader_writer_test
    LIBS parallel_no_concurrency dca_io function)
  if (DCA_HAVE_MPI)
    dca_add_gtest(adios2_parallelIO_test
      MPI MPI_NUMPROC 12
      LIBS dca_io parallel_mpi_concurrency function parallel_util)
  endif()
endif()
