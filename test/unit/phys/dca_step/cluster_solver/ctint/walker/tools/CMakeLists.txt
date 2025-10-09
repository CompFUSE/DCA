set(TEST_INCLUDES ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
set(TEST_LIBS     ${DCA_LIBS})

set(THIS_TEST_LIBS ${DCA_LIBS} ctint)

dca_add_gtest(d_matrix_builder_gpu_test
    FAST
    CUDA
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${THIS_TEST_LIBS} g0_interpolation ctint ${DCA_KERNEL_LIBS}
    )
