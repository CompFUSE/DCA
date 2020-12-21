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

dca_add_gtest(matrix_view_test GTEST_MAIN)

dca_add_gtest(matrix_cpu_gpu_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS})

dca_add_gtest(matrix_gpu_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS})

dca_add_gtest(matrixop_complex_cpu_test
              GTEST_MAIN
              LIBS ${DCA_LIBS})

dca_add_gtest(matrixop_real_cpu_test
              GTEST_MAIN
              LIBS ${DCA_LIBS})

dca_add_gtest(matrixop_cpu_gpu_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS})

dca_add_gtest(matrixop_real_gpu_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS})

dca_add_gtest(matrixop_complex_gpu_test
              GTEST_MAIN
              CUDA
              LIBS ${DCA_LIBS})

dca_add_gtest(reshapable_matrix_cpu_test GTEST_MAIN)

dca_add_gtest(reshapable_matrix_gpu_test CUDA GTEST_MAIN LIBS cuda_utils)

dca_add_gtest(multi_vector_test GTEST_MAIN)

add_subdirectory(util)
