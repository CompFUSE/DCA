set(TEST_INCLUDES ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
set(TEST_LIBS     ${DCA_LIBS})

dca_add_gtest(ct_int_interpolation_test
    FAST
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )

# disable due to function.hpp inclusion in cuda code.
#
# dca_add_gtest(ct_int_interpolation_gpu_test
#     FAST
#     CUDA
#     GTEST_MAIN
#     INCLUDE_DIRS ${TEST_INCLUDES}
#     LIBS     ${TEST_LIBS}
#     )

set(THIS_TEST_LIBS ${DCA_LIBS} mc_kernels; ctint)
if (DCA_HAVE_ADIOS2)
  set(THIS_TEST_LIBS ${THIS_TEST_LIBS} dca_adios2)
endif()


dca_add_gtest(d_matrix_builder_gpu_test
    FAST
    CUDA
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${THIS_TEST_LIBS}
    )
