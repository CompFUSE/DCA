# linalg unit tests

dca_add_gtest(vector_cpu_test GTEST_MAIN)

dca_add_gtest(vector_cpu_gpu_test
              GTEST_MAIN
              CUDA
              INCLUDE_DIRS ${TEST_INCLUDE_DIRS}
              LIBS ${DCA_LIBS})

dca_add_gtest(vector_gpu_test
              GTEST_MAIN
              CUDA
              INCLUDE_DIRS ${TEST_INCLUDE_DIRS}
              LIBS ${DCA_LIBS})

dca_add_gtest(matrix_cpu_test GTEST_MAIN)

dca_add_gtest(matrix_cpu_gpu_test
              GTEST_MAIN
              CUDA
              INCLUDE_DIRS ${TEST_INCLUDE_DIRS}
              LIBS ${DCA_LIBS})

dca_add_gtest(matrix_gpu_test
              GTEST_MAIN
              CUDA
              INCLUDE_DIRS ${TEST_INCLUDE_DIRS}
              LIBS ${DCA_LIBS})
