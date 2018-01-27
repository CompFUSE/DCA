# Pthreading unit tests

dca_add_gtest(threading_test
  GTEST_MAIN
  STDTHREAD
  LIBS ${DCA_THREADING_LIB})
