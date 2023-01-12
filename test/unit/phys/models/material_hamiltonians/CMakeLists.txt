# Material Hamiltonians unit tests

dca_add_gtest(material_lattice_NiO_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function ${LAPACK_LIBRARIES} ${DCA_GPU_LIBS} ${UTIL_LIBS} models)

dca_add_gtest(material_lattice_FeSn_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS function ${LAPACK_LIBRARIES} ${DCA_GPU_LIBS} ${UTIL_LIBS} models)

set_tests_properties(material_lattice_FeSn_test PROPERTIES WILL_FAIL 1)
