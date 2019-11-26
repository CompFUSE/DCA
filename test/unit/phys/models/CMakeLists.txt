add_subdirectory(analytic_hamiltonians)
add_subdirectory(material_hamiltonians)

dca_add_gtest(general_interaction_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIRS}
  LIBS function ${LAPACK_LIBRARIES} ${DCA_CUDA_LIBS})
