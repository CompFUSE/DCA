
find_library(PAPI_LIB papi)
if(PAPI_LIB AND DCA_PROFILER EQUAL "PAPI")
    dca_add_gtest(papi_profiler_test GTEST_MAIN LIBS json profiling)
endif()
