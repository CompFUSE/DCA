# test/unit/phys/dca_step/cluster_solver/shared_tools/sp

dca_add_gtest(sp_accumulator_gpu_test
  CUDA
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     FFTW::Double ${DCA_LIBS}
  )

dca_add_gtest(sp_accumulator_single_meas_G_test
  CUDA
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     FFTW::Double ${DCA_LIBS}
  )

#;mc_kernels
