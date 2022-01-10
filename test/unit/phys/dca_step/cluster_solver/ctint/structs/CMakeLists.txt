set(TEST_INCLUDES ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
#set(TEST_LIBS     ${DCA_LIBS})

dca_add_gtest(ct_int_interaction_vertices_test
    FAST
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    CUDA
    LIBS ${DCA_LIBS}
    )

dca_add_gtest(solver_configuration_test
    FAST
    GTEST_MAIN
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS random;ctint;time_and_frequency_domains;function;json;${DCA_THREADING_LIBS};dca_io
    )

dca_add_gtest(read_write_ctint_config_test
  FAST
  GTEST_MAIN
  INCLUDE_DIRS ${TEST_INCLUDES}
  LIBS random;time_and_frequency_domains;function;ctint;${DCA_THREADING_LIBS};dca_io
  )

if (DCA_HAVE_GPU)
  target_link_libraries(read_write_ctint_config_test  PRIVATE gpu_utils)
endif()
