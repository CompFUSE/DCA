# test/unit/phys/dca_loop/disorder

# DISORDERED_G0 numeric kernels: translationalAverage orientation round-trip (Nc=1, 2x2, 3x4)
# and the per-frequency disorder Dyson (accumulateDisorderDyson). Uses mock domains, so it needs
# no FFTW and no concurrency-initialized cluster domain.
dca_add_gtest(disorder_average_test
  FAST
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS function ${DCA_LIBS} parallel_no_concurrency)

# DISORDERED_G0 disorder-averaging pipeline (loop body of workTheClusters + averageGkw) on the real
# DcaData with U=0 (M=0): clean-limit identity and the weighted disorder average. Uses the real
# cluster (2x4) and the r<->k FunctionTransform, so it links FFTW.
dca_add_gtest(disorder_loop_test
  FAST
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIR}
  LIBS function FFTW::Double ${DCA_LIBS} parallel_no_concurrency)
