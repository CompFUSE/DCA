# stdthread unit tests

dca_add_gtest(stdthread_test GTEST_MAIN LIBS ${DCA_THREADING_LIBS})

add_subdirectory(thread_pool)
