set(TEST_INCLUDES ${DCA_INCLUDE_DIRS})
set(TEST_LIBS     ${DCA_LIBS})

dca_add_gtest(stdthread_ctaux_gpu_tp_cluster_solver_test
    EXTENSIVE
    CUDA
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )

