# DCA algorithms unit tests

set(gtest_main_flag GTEST_MAIN)
if(DCA_HAVE_ADIOS2)
  set(gtest_main_flag GTEST_MPI_MAIN)
endif(DCA_HAVE_ADIOS2)

dca_add_gtest(compute_band_structure_test
  GTEST_MPI_MAIN
  INCLUDE_DIRS ${PROJECT_SOURCE_DIR}
  LIBS FFTW::Double dca_io json function cluster_domains enumerations time_and_frequency_domains quantum_domains ${UTIL_LIBS}
  ${LAPACK_LIBRARIES} parallel_no_concurrency)

if (DCA_HAVE_ADIOS2)
  set(optional_io_lib dca_io)
endif()

# \todo I think this should be able to run without GPU
dca_add_gtest(compute_free_greens_function_test
  ${gtest_main_flag}
  CUDA
  LIBS function FFTW::Double ${DCA_THREADING_LIBS} ${LAPACK_LIBRARIES} ${UTIL_LIBS} ${DCA_GPU_LIBS} ${optional_io_lib} parallel_util parallel_no_concurrency)

dca_add_gtest(compute_greens_function_test
  ${gtest_main_flag}
  LIBS function FFTW::Double ${DCA_THREADING_LIBS} ${LAPACK_LIBRARIES} ${UTIL_LIBS} ${DCA_GPU_LIBS} ${optional_io_lib} parallel_util parallel_no_concurrency)
