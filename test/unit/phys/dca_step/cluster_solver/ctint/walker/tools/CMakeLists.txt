set(TEST_INCLUDES ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
set(TEST_LIBS     ${DCA_LIBS})

set(THIS_TEST_LIBS ${DCA_LIBS} ctint)
if (DCA_HAVE_ADIOS2)
  set(THIS_TEST_LIBS ${THIS_TEST_LIBS} dca_adios2)
endif()

dca_add_gtest(d_matrix_builder_gpu_test
    FAST
    CUDA
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${THIS_TEST_LIBS} g0_interpolation
    )
