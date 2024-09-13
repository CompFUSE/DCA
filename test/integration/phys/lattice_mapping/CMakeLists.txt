# Lattice mapping tests

dca_add_gtest(lattice_mapping_sp_test EXTENSIVE
  GTEST_MAIN
  CUDA
  INCLUDE_DIRS ${SIMPLEX_GM_RULE_INCLUDE_DIR};${FFTW_INCLUDE_DIR};${PROJECT_SOURCE_DIR}
  LIBS json function function_transform cluster_domains ${DCA_THREADING_LIBS} dca_io parallel_util
       time_and_frequency_domains
       quantum_domains gaussian_quadrature
       tetrahedron_mesh coarsegraining enumerations dca_hdf5 ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES} ${DCA_GPU_LIBS}
       lapack gpu_utils parallel_no_concurrency)
