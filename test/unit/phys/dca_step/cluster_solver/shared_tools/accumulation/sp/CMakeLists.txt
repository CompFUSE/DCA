# test/unit/phys/dca_step/cluster_solver/shared_tools/sp

dca_add_gtest(sp_accumulator_gpu_test
  CUDA
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS}
  )
