# nfft tests

set(THIS_TEST_LIBS time_and_frequency_domains random function mc_kernels gpu_utils nfft)
if (DCA_HAVE_ADIOS2)
  set(THIS_TEST_LIBS ${THIS_TEST_LIBS} dca_adios2)
else()
  set(THIS_TEST_LIBS ${THIS_TEST_LIBS} hdf5::hdf5_cpp)
endif()

message("thread option for nfft test: ${test_thread_option}")

include(CMakePrintHelpers)
cmake_print_properties(TARGETS HPX::hpx
  PROPERTIES IMPORTED INTERFACE_LINK_LIBRARIES IMPORTED_LINK_FLAGS LOCATION)

dca_add_gtest(dnfft_1d_gpu_test
  THREADED
  GTEST_MAIN
  CUDA
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR};${PROJECT_SOURCE_DIR}
  LIBS ${FFTW_LIBRARY} time_and_frequency_domains random function mc_kernels gpu_utils nfft)

cmake_print_properties(TARGETS dnfft_1d_gpu_test
  PROPERTIES LINK_LIBRARIES)
