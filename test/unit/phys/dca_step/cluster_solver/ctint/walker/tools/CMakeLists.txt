dca_add_gtest(d_matrix_builder_gpu_test
    FAST
    CUDA
    GTEST_MAIN
    INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
    LIBS     ${DCA_LIBS}
    )
