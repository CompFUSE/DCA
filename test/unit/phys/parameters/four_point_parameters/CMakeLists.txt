# Four-point parameters' unit tests

set(FOUR_POINT_PARAMETERS_LIBS json enumerations ${DCA_GPU_LIBS})

# Conditional linking for MAGMA (only link if MAGMA is available)
if (MAGMA_FOUND)
# This is only required for cuda 12.8 but harmless for earlier cuda
# and rocm
# why this test pulls in a bunch of magma references is a mystery
  list(APPEND FOUR_POINT_PARAMETERS_LIBS magma::magma)
endif()

# Create the test, passing the library list via the LIBS option
dca_add_gtest(four_point_parameters_test
  GTEST_MAIN
  LIBS ${FOUR_POINT_PARAMETERS_LIBS})
