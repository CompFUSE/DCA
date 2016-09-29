add_subdirectory(domains)

dca_add_gtest(function_library_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
