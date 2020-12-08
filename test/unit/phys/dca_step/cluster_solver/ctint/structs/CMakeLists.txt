set(TEST_INCLUDES ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
#set(TEST_LIBS     ${DCA_LIBS})

dca_add_gtest(ct_int_interaction_vertices_test
    FAST
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS ${DCA_LIBS}
    )

dca_add_gtest(solver_configuration_test
    FAST
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS random;ctint;time_and_frequency_domains;function;json
    )

dca_add_gtest(read_write_ctint_config_test
  FAST
  GTEST_MAIN
  INCLUDE_DIRS ${TEST_INCLUDES}
  LIBS random;time_and_frequency_domains;function;ctint
  )
