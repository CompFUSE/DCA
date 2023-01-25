# Model parameters' unit tests

# dca_add_gtest(model_parameters_2d_2band_test
#   GTEST_MAIN
#   LIBS json)

# dca_add_gtest(model_parameters_2d_4band_test
#   GTEST_MAIN
#   LIBS json)

dca_add_gtest(model_parameters_bilayer_hubbard_test
  GTEST_MAIN
  LIBS function_transform json)

dca_add_gtest(model_parameters_material_test
  GTEST_MAIN
  LIBS function_transform json models)

dca_add_gtest(model_parameters_single_band_hubbard_test
  GTEST_MAIN
  LIBS function_transform json)
