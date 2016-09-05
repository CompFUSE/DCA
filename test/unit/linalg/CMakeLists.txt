# linalg unit tests

add_gtest(NAME vector_cpu_test
          GTEST_MAIN)

add_gtest(NAME vector_cpu_gpu_test
          CUDA
          INCLUDES ${DCA_INCLUDES}
          LIBS ${DCA_LIBRARIES}
          GTEST_MAIN)

add_gtest(NAME vector_gpu_test
          CUDA
          INCLUDES ${DCA_INCLUDES}
          LIBS ${DCA_LIBRARIES}
          GTEST_MAIN)

add_gtest(NAME matrix_cpu_test
          GTEST_MAIN)

add_gtest(NAME matrix_cpu_gpu_test
          CUDA
          INCLUDES ${DCA_INCLUDES}
          LIBS ${DCA_LIBRARIES}
          GTEST_MAIN)

add_gtest(NAME matrix_gpu_test
          CUDA
          INCLUDES ${DCA_INCLUDES}
          LIBS ${DCA_LIBRARIES}
          GTEST_MAIN)
