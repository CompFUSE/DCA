# statistical_testing integration test

dca_add_gtest(statistical_testing_integration_test
    GTEST_MAIN
    EXTENSIVE
    INCLUDE_DIRS
    LIBS ${LAPACK_LIBRARIES};statistical_testing)

