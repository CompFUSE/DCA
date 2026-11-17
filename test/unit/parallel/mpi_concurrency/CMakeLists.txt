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
# Exercise the full suite at several sizes without compiling additional executables.
# Keep the original test name for the six-rank run.
if(TARGET mpi_collective_sum_test)
  set_tests_properties(mpi_collective_sum_test PROPERTIES PROCESSORS 6)
  foreach(numproc IN ITEMS 1 2 3)
    if(numproc EQUAL 1)
      set(test_suffix _1_rank)
    else()
      set(test_suffix _${numproc}_ranks)
    endif()
    add_test(NAME mpi_collective_sum_test${test_suffix}
      COMMAND ${TEST_RUNNER} ${MPIEXEC_NUMPROC_FLAG} ${numproc}
              ${MPIEXEC_PREFLAGS} ${SMPIARGS_FLAG_MPI} ${MPIEXEC_POSTFLAGS} ${CVD_LAUNCHER}
              "$<TARGET_FILE:mpi_collective_sum_test>" ${DCA_TESTING_FLAGS})
    set_tests_properties(mpi_collective_sum_test${test_suffix} PROPERTIES PROCESSORS ${numproc})
  endforeach()
endif()
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
