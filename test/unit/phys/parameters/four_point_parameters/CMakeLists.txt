# Four-point parameters' unit tests

# why this test pulls in a bunch of magma references is a mystery
dca_add_gtest(four_point_parameters_test
  GTEST_MAIN
  LIBS json enumerations ${DCA_GPU_LIBS})
# Conditional linking for MAGMA (only link if MAGMA is available)
if (MAGMA_FOUND)
  target_link_libraries(four_point_parameters_test PRIVATE magma::magma)
endif()
