add_subdirectory(analytic_hamiltonians)
add_subdirectory(material_hamiltonians)

dca_add_gtest(general_interaction_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function function_transform ${LAPACK_LIBRARIES} ${DCA_GPU_LIBS} ${UTIL_LIBS})

dca_add_gtest(traits_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function function_transform ${LAPACK_LIBRARIES} ${DCA_GPU_LIBS} ${UTIL_LIBS})
