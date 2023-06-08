# IO's unit tests

dca_add_gtest(dca_data_test
  FAST
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS function ${DCA_LIBS})
