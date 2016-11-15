# Model parameters' unit tests

dca_add_gtest(model_parameters_2d_2band_test
  GTEST_MAIN
  LIBS json)

dca_add_gtest(model_parameters_2d_4band_test
  GTEST_MAIN
  LIBS json)

dca_add_gtest(model_parameters_2d_bilayer_test
  GTEST_MAIN
  LIBS json)

dca_add_gtest(model_parameters_material_test
  GTEST_MAIN
  LIBS json)

dca_add_gtest(model_parameters_tight_binding_test
  GTEST_MAIN
  LIBS json)
