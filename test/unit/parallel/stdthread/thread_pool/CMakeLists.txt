# thread pool unit tests

dca_add_gtest(thread_pool_test GTEST_MAIN LIBS parallel_stdthread)

EXECUTE_PROCESS( COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE )

# temporary fix, we need some platform and compiler awareness at a highler level in our cmake
# this test also depends on GNU compiler extensions, something like hwloc should be used.
IF(NOT (${ARCHITECTURE} STREQUAL "ppc64le"))
  dca_add_gtest(affinity_test GTEST_MAIN LIBS parallel_stdthread)
ENDIF()
