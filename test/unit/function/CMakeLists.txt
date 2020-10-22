# Function library's unit tests

add_subdirectory(domains)
add_subdirectory(util)

dca_add_gtest(function_library_test
  GTEST_MAIN
  LIBS function)

if (function_library_test)
  # We want to test self-assignment and self-move.
  target_compile_options(function_library_test PRIVATE -Wno-self-assign-overloaded -Wno-self-move)
endif()

dca_add_gtest(set_to_zero_test GTEST_MAIN)

dca_add_gtest(function_construction_test
  GTEST_MAIN
  MPI MPI_NUMPROC 4
  LIBS parallel_mpi_concurrency;function)
