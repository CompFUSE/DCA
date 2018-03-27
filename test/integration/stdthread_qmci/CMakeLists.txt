set(TEST_INCLUDES ${DCA_INCLUDE_DIRS})
set(TEST_LIBS     ${DCA_LIBS})

dca_add_gtest(stdthread_ctaux_tp_cluster_solver_test
    EXTENSIVE
    CUDA
    GTEST_MAIN
    INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
    LIBS ${DCA_LIBS})

