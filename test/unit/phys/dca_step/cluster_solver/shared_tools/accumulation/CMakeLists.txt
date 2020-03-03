dca_add_gtest(time_correlator_gpu_test
        CUDA
        GTEST_MAIN
        INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
        LIBS ${DCA_LIBS})

dca_add_gtest(time_correlator_cpu_test
        GTEST_MAIN
        INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
        LIBS ${DCA_LIBS})

add_subdirectory(sp)
add_subdirectory(tp)
