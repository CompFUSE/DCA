set(TEST_INCLUDES ${DCA_INCLUDE_DIRS})
set(TEST_LIBS     ${DCA_LIBS})

dca_add_gtest(posix_ctaux_cluster_solver_test
    EXTENSIVE
    GTEST_MAIN
    PTHREADS
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )

