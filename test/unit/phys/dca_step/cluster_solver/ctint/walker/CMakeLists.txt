set(TEST_INCLUDES ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
set(TEST_LIBS     ${DCA_LIBS})

dca_add_gtest(ct_int_interpolation_test
    FAST
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )

dca_add_gtest(ct_int_interpolation_test_gpu
    FAST
    CUDA
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS};ctint_kernels
    )

dca_add_gtest(ct_int_walker_test
    FAST
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )

dca_add_gtest(ct_int_ndd_walker_test
    FAST
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )

dca_add_gtest(ct_int_walker_submatrix_test
        FAST
        GTEST_MAIN
        INCLUDE_DIRS ${TEST_INCLUDES}
        LIBS     ${TEST_LIBS}
        )

dca_add_gtest(ct_int_walker_submatrix_gpu_test
        CUDA
        FAST
        GTEST_MAIN
        INCLUDE_DIRS ${TEST_INCLUDES}
        LIBS     ${TEST_LIBS}
        )

dca_add_gtest(d_matrix_builder_test_gpu
    FAST
    CUDA
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS};ctint_kernels;ctint
    )
