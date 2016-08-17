# MCI parameters unit tests

set(TEST_INCLUDES ${PROJECT_SOURCE_DIR}/src;${FFTW_INCLUDE_DIR})
set(TEST_LIBS ${LAPACK_LIBRARIES})

add_gtest(NAME mci_parameters_test
  GTEST_MAIN
  INCLUDES ${TEST_INCLUDES}
  LIBS ${TEST_LIBS})
target_compile_definitions(mci_parameters_test PRIVATE DCA_TEST_DIRECTORY="${CMAKE_CURRENT_SOURCE_DIR}/")
