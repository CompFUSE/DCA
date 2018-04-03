# Inference unit tests

dca_add_gtest(richardson_lucy_deconvolution_test
  GTEST_MAIN
  LIBS function ${LAPACK_LIBRARIES})
