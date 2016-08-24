set(TEST_INCLUDES "${PROJECT_SOURCE_DIR}/include;${PROJECT_SOURCE_DIR}/src;${FFTW_INCLUDE_DIR}")

add_gtest(NAME collective_sum_interface_mpi_test
          MPI
          MPI_NUMPROC 8
          INCLUDES "${TEST_INCLUDES}"
          LIBS ${FFTW_LIBRARY})
