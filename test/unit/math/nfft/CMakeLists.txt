# nfft tests

set(THIS_TEST_LIBS time_and_frequency_domains random function mc_kernels cuda_utils nfft)
if (DCA_HAVE_ADIOS2)
  set(THIS_TEST_LIBS ${THIS_TEST_LIBS} dca_adios2)
endif()

dca_add_gtest(dnfft_1d_gpu_test
  HPX
  GTEST_MAIN
  CUDA
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR};${PROJECT_SOURCE_DIR}
  LIBS ${FFTW_LIBRARY} time_and_frequency_domains random function mc_kernels cuda_utils nfft)
