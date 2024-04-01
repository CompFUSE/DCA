# (Phys) deconvolution unit tests

dca_add_gtest(deconvolution_routines_test
  GTEST_MPI_MAIN
  CUDA
  INCLUDE_DIRS ${SIMPLEX_GM_RULE_INCLUDE_DIR} ${FFTW_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}
  LIBS dca_io function function_transform cluster_domains  ${DCA_THREADING_LIBS} parallel_util parallel_no_concurrency time_and_frequency_domains
       quantum_domains gaussian_quadrature tetrahedron_mesh coarsegraining enumerations
       ${LAPACK_LIBRARIES} lapack ${UTIL_LIBS} )
