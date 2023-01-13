# nfft tests

set(THIS_TEST_LIBS time_and_frequency_domains random function mc_kernels gpu_utils nfft)
if (DCA_HAVE_ADIOS2)
  set(THIS_TEST_LIBS ${THIS_TEST_LIBS} dca_adios2)
else()
  set(THIS_TEST_LIBS ${THIS_TEST_LIBS} hdf5::hdf5_cpp)
endif()

dca_add_gtest(dnfft_1d_gpu_test
  THREADED
  GTEST_MAIN
  CUDA
  INCLUDE_DIRS ${PROJECT_SOURCE_DIR}
  LIBS FFTW::Double time_and_frequency_domains random function mc_kernels gpu_utils nfft)
