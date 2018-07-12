# stdthread unit tests

dca_add_gtest(stdthread_test GTEST_MAIN LIBS parallel_stdthread)

add_subdirectory(thread_pool)
