# Four-point parameters' unit tests

dca_add_gtest(four_point_parameters_test
  GTEST_MAIN
  LIBS json enumerations ${DCA_CUDA_LIBS})
