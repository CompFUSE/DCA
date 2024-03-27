set(TEST_INCLUDES ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
set(TEST_LIBS     ${DCA_LIBS})

if(DCA_HAVE_GPU)

message("DCA_KERNEL_LIBS: ${DCA_KERNEL_LIBS}")
  
dca_add_gtest(ct_aux_walker_test
    FAST
    GTEST_MPI_MAIN
    CUDA
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     FFTW::Double ${TEST_LIBS} g0_interpolation ${DCA_KERNEL_LIBS}
    )
endif()
