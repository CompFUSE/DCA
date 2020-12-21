set(TEST_INCLUDES ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
set(TEST_LIBS     ${DCA_LIBS})

dca_add_gtest(ct_int_interpolation_test
    FAST
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )

dca_add_gtest(ct_int_interpolation_gpu_test
    FAST
    CUDA
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )

dca_add_gtest(d_matrix_builder_gpu_test
    FAST
    CUDA
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )
