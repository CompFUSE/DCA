# test/unit/phys/dca_step/cluster_solver/shared_tools

dca_add_gtest(tp_accumulator_test
  FAST
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS ${DCA_LIBS})

add_subdirectory(tp/ndft)
