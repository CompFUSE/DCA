# nfft tests

dca_add_gtest(dnfft_1d_gpu_test
  GTEST_MAIN
  CUDA
  INCLUDE_DIRS ${FFTW_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS ${FFTW_LIBRARIES} time_and_frequency_domains random function dnfft_kernels cuda_utils nfft)
