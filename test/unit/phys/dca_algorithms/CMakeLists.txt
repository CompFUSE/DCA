# DCA algorithms unit tests

dca_add_gtest(compute_band_structure_test
  GTEST_MAIN
  INCLUDE_DIRS ${PROJECT_SOURCE_DIR}
  LIBS FFTW::Double dca_io json function cluster_domains enumerations time_and_frequency_domains quantum_domains ${UTIL_LIBS}
  ${LAPACK_LIBRARIES})

if (DCA_HAVE_ADIOS2)
  set(optional_io_lib dca_io)
endif()

# \todo I think this should be able to run without GPU
dca_add_gtest(compute_free_greens_function_test
  GTEST_MAIN
  CUDA
  LIBS function FFTW::Double ${DCA_THREADING_LIBS} ${LAPACK_LIBRARIES} ${UTIL_LIBS} ${DCA_GPU_LIBS} ${optional_io_lib} parallel_util)

dca_add_gtest(compute_greens_function_test
  GTEST_MAIN
  LIBS function FFTW::Double ${DCA_THREADING_LIBS} ${LAPACK_LIBRARIES} ${UTIL_LIBS} ${DCA_GPU_LIBS} ${optional_io_lib} parallel_util)
