# unit/math/function_transform

set(gtest_main_flag GTEST_MAIN)
if(DCA_HAVE_ADIOS2)
  set(gtest_main_flag GTEST_MPI_MAIN)
endif(DCA_HAVE_ADIOS2)

dca_add_gtest(function_transform_test
        ${gtest_main_flag}
        INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIR}
        LIBS FFTW::Double ${DCA_LIBS} parallel_no_concurrency)

dca_add_gtest(space_transform_2D_test
        INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIR}
        LIBS FFTW::Double ${DCA_LIBS} parallel_no_concurrency)

dca_add_gtest(space_transform_2D_gpu_test
    CUDA
    INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${PROJECT_SOURCE_DIR}
    LIBS dca_io FFTW::Double ${DCA_LIBS})# magma::magma magma::sparse)
