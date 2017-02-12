add_subdirectory(analytic_hamiltonians)
add_subdirectory(material_hamiltonians)

dca_add_gtest(general_interaction_test
  GTEST_MAIN
  LIBS function ${LAPACK_LIBRARIES})
