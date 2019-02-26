# test/unit/phys/dca_step/cluster_solver/shared_tools/tp

dca_add_gtest(tp_accumulator_test
  FAST
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS}
  )

dca_add_gtest(tp_accumulator_singleband_test
  FAST
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS}
  )

dca_add_gtest(tp_accumulator_gpu_test
  CUDA
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS}
  )

dca_add_gtest(tp_accumulator_singleband_gpu_test
  CUDA
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS}
  )

dca_add_gtest(g4_name_to_four_point_type_test GTEST_MAIN)

add_subdirectory(ndft)
