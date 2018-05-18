dca_add_gtest(function_cut_test
  GTEST_MAIN
    INCLUDES ${DCA_INCLUDE_DIRS}
    LIBS ${DCA_LIBS})

dca_add_gtest(statistical_testing_test
    GTEST_MAIN
    INCLUDES ${DCA_INCLUDE_DIRS}
    LIBS ${LAPACK_LIBRARIES} ${DCA_LIBS} statistical_testing)
