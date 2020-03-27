set(TEST_INCLUDES ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
set(TEST_LIBS     ${DCA_LIBS})

dca_add_gtest(ctint_square_lattice_test
    EXTENSIVE
    MPI MPI_NUMPROC 4
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )

dca_add_gtest(ctint_square_lattice_test_gpu
    EXTENSIVE
    CUDA
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS};mc_kernels;ctint;mc_kernels
    )

dca_add_gtest(ctint_hund_lattice_test
    EXTENSIVE
    MPI MPI_NUMPROC 4
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )

dca_add_gtest(ctint_fe_as_test
    EXTENSIVE
    MPI MPI_NUMPROC 4
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )

dca_add_gtest(ctint_double_update_comparison_test
        EXTENSIVE
        GTEST_MAIN
        INCLUDE_DIRS ${TEST_INCLUDES}
        LIBS     ${TEST_LIBS}
   )

dca_add_gtest(ctint_square_lattice_tp_test
    EXTENSIVE
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )
