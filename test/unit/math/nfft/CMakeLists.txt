# nfft tests

dca_add_gtest(dnfft_1d_gpu_test
  HPX
  GTEST_MAIN
  CUDA
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR};${PROJECT_SOURCE_DIR}
  LIBS ${FFTW_LIBRARY} time_and_frequency_domains random function dnfft_kernels cuda_utils nfft ${DCA_LIBS})
