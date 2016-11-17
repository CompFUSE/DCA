# DCA algrotihms unit tests

dca_add_gtest(compute_band_structure_test
  GTEST_MAIN
  INCLUDE_DIRS ${FFTW_INCLUDE_DIR}
  LIBS dca_algorithms json_reader cluster_domains quantum_domains time_and_frequency_domains
  function ${FFTW_LIBRARY})
