# Tetrahedron mesh unit tests

dca_add_gtest(facet_test
  GTEST_MAIN
  LIBS tetrahedron_mesh)

dca_add_gtest(tetrahedron_test
  GTEST_MAIN
  LIBS tetrahedron_mesh)
