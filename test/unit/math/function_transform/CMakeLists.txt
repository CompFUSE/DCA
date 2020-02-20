# unit/math/function_transform
dca_add_gtest(function_transform_test
        GTEST_MAIN
        INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIRS}
        LIBS ${DCA_LIBS})

dca_add_gtest(space_transform_2D_test
        GTEST_MAIN
        INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIRS}
        LIBS ${DCA_LIBS})


dca_add_gtest(space_transform_2D_gpu_test
    GTEST_MAIN
    CUDA
    INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIRS}
    LIBS ${DCA_LIBS})
