dca_add_gtest(function_cut_test
    GTEST_MAIN
    INCLUDES ${DCA_INCLUDES}
    LIBS ${DCA_LIBS})

dca_add_gtest(statistical_testing_test
    GTEST_MAIN
    INCLUDES ${DCA_INCLUDES}
    LIBS ${LAPACK_LIBRARIES};statistical_testing)
