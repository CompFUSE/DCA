# test/unit/phys/dca_step/cluster_solver/shared_tools/tp/ndft

# dca_add_gtest(cached_ndft_cpu_test
#   GTEST_MAIN
#   CUDA
#   FAST
#   INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIRS}
#   LIBS ${LAPACK_LIBRARIES} time_and_frequency_domains random function dca_io ${FFTW_LIBRARIES})

set(THIS_TEST_LIBS time_and_frequency_domains random function gpu_utils)
if (DCA_HAVE_ADIOS2)
  set(THIS_TEST_LIBS ${THIS_TEST_LIBS} dca_adios2)
endif()

dca_add_gtest(cached_ndft_gpu_test
  THREADED
  GTEST_MAIN
  CUDA
  FAST
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS function_transform ${DCA_KERNEL_LIBS} hdf5::hdf5_cpp magma::magma magma::sparse ${THIS_TEST_LIBS})
