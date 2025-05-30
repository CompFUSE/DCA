add_subdirectory(tools)

set(TEST_INCLUDES ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
set(TEST_LIBS     ${DCA_LIBS})

dca_add_gtest(ct_int_walker_test
    FAST
    GTEST_MPI_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     FFTW::Double ${TEST_LIBS} g0_interpolation ${DCA_KERNEL_LIBS}
    )

dca_add_gtest(ct_int_walker_submatrix_test
        FAST
        GTEST_MPI_MAIN
        INCLUDE_DIRS ${TEST_INCLUDES}
        LIBS     ${TEST_LIBS} g0_interpolation ${DCA_GPU_LIBS}
        )

dca_add_gtest(ct_int_walker_multinsertion_test
        FAST
        GTEST_MPI_MAIN
        INCLUDE_DIRS ${TEST_INCLUDES}
        LIBS     ${TEST_LIBS} g0_interpolation
        )

dca_add_gtest(ct_int_walker_submatrix_gpu_test
        CUDA
        FAST
        GTEST_MPI_MAIN
        INCLUDE_DIRS ${TEST_INCLUDES}
        LIBS     ${TEST_LIBS} g0_interpolation ${DCA_GPU_LIBS} ctint ${DCA_KERNEL_LIBS}
      )

dca_add_gtest(ct_int_walker_submatrix_gpu_complex_test
        CUDA
        FAST
        GTEST_MPI_MAIN
        INCLUDE_DIRS ${TEST_INCLUDES}
        LIBS     ${TEST_LIBS} g0_interpolation ${DCA_GPU_LIBS} ctint ${DCA_KERNEL_LIBS}
        )
