# unit/math/function_transform
dca_add_gtest(function_transform_test
        INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIR}
        LIBS FFTW::Double ${DCA_LIBS})

dca_add_gtest(space_transform_2D_test
        INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIR}
        LIBS FFTW::Double ${DCA_LIBS})


dca_add_gtest(space_transform_2D_gpu_test
    GTEST_MAIN
    CUDA
    INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${PROJECT_SOURCE_DIR}
    LIBS FFTW::Double ${DCA_LIBS})# magma::magma magma::sparse)
