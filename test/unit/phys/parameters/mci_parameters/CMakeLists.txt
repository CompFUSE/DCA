# MCI parameters unit tests

set(TEST_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/src ${FFTW_INCLUDE_DIR})
set(TEST_LIBS ${LAPACK_LIBRARIES})

dca_add_gtest(mci_parameters_test
  GTEST_MAIN
  INCLUDE_DIRS ${TEST_INCLUDE_DIRS}
  LIBS ${TEST_LIBS})
