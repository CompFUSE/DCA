# Analytic Hamiltonians unit tests

dca_add_gtest(bilayer_lattice_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function ${LAPACK_LIBRARIES} ${DCA_CUDA_LIBS})

dca_add_gtest(square_lattice_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function ${LAPACK_LIBRARIES} ${DCA_CUDA_LIBS})

dca_add_gtest(triangular_lattice_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function ${LAPACK_LIBRARIES} ${DCA_CUDA_LIBS})
