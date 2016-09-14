# linalg unit tests

dca_add_gtest(vector_cpu_test GTEST_MAIN)

dca_add_gtest(vector_cpu_gpu_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS})

dca_add_gtest(vector_gpu_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS})

dca_add_gtest(matrix_cpu_test GTEST_MAIN)

dca_add_gtest(matrix_cpu_gpu_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS})

dca_add_gtest(matrix_gpu_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS})

dca_add_gtest(matrixop_cpu_test
              GTEST_MAIN
              LIBS ${DCA_LIBS})

dca_add_gtest(matrixop_gpu_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS})
