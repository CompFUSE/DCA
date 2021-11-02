# thread pool unit tests
if(DCA_HAVE_HPX)
dca_add_gtest(hpx_thread_pool_test GTEST_MAIN HPX LIBS ${DCA_THREAD_LIBS})
#dca_add_gtest(affinity_test GTEST_MAIN LIBS ${DCA_THREADING_LIBS})
endif()
