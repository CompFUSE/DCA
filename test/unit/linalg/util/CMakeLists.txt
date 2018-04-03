# linalg/util unit tests

dca_add_gtest(check_RC_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS} ${DCA_CUDA_LIBS})

dca_add_gtest(stream_handle_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS})
