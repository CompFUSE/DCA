dca_add_gtest(g0_interpolation_test
        FAST
        GTEST_MAIN
        INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
        LIBS  g0_interpolation
        )

dca_add_gtest(g0_interpolation_gpu_test
        FAST
        CUDA
        GTEST_MAIN
        INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
        LIBS  mc_kernels g0_interpolation
        )
