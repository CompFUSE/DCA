# statistical_testing integration test

dca_add_gtest(statistical_testing_integration_test
    GTEST_MAIN
    EXTENSIVE
    INCLUDE_DIRS
    LIBS ${LAPACK_LIBRARIES} ${DCA_CUDA_LIBS} statistical_testing)

