# statistical_testing integration test

message("statistical_testing integration test libs: ${LAPACK_LIBRARIES} ${DCA_GPU_LIBS} statistical_testing gpu_utils")

set(UTIL_LIBS)
if(DCA_HAVE_GPU)
  set(UTIL_LIBS gpu_utils)
endif()

dca_add_gtest(statistical_testing_integration_test
    GTEST_MAIN
    EXTENSIVE
    INCLUDE_DIRS
    LIBS ${LAPACK_LIBRARIES} ${DCA_GPU_LIBS} statistical_testing ${UTIL_LIBS})

