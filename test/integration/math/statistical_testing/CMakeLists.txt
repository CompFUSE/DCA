# statistical_testing integration test

message("statistical_testing integration test libs: ${LAPACK_LIBRARIES} ${DCA_GPU_LIBS} statistical_testing gpu_utils")

dca_add_gtest(statistical_testing_integration_test
    GTEST_MAIN
    EXTENSIVE
    INCLUDE_DIRS
    LIBS ${LAPACK_LIBRARIES} ${DCA_GPU_LIBS} statistical_testing gpu_utils)

