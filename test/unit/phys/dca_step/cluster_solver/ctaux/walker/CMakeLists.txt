set(TEST_INCLUDES ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
set(TEST_LIBS     ${DCA_LIBS})

if(DCA_HAVE_GPU)
dca_add_gtest(ct_aux_walker_test
    FAST
    GTEST_MAIN
    CUDA
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     FFTW::Double ${TEST_LIBS} g0_interpolation ${KERNELS_LIB}
    )
endif()
