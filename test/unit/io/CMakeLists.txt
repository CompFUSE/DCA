# IO's unit tests

dca_add_gtest(buffer_test
  GTEST_MAIN)

dca_add_gtest(reader_test
  GTEST_MAIN
  LIBS function dca_hdf5 json hdf5::hdf5_cpp)

dca_add_gtest(writer_test
  GTEST_MAIN
  LIBS function dca_hdf5 json hdf5::hdf5_cpp)

dca_add_gtest(json_reader_test
        GTEST_MAIN
        LIBS json)

dca_add_gtest(reader_writer_test
  LIBS dca_hdf5 hdf5::hdf5_cpp json function)

if (DCA_HAVE_ADIOS2)
  dca_add_gtest(adios2_reader_writer_test
    GTEST_MAIN
    LIBS parallel_no_concurrency dca_adios2 adios2::adios2 function)
  if (DCA_HAVE_MPI)
    dca_add_gtest(adios2_parallelIO_test
      MPI MPI_NUMPROC 12
      LIBS dca_adios2 adios2::adios2 parallel_mpi_concurrency;function;parallel_util)
  endif()
endif()
