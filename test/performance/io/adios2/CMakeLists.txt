# Adios2 parallel write performance tests

include(CMakeParseArguments)

if (DCA_WITH_TESTS_PERFORMANCE)

function(adios_scaling_test case num_ranks)

  set(CASE_NAME "${case}_${num_ranks}")
  add_executable(${CASE_NAME} adios2_scalable_test.cpp)
  target_link_libraries(${CASE_NAME} PRIVATE ${DCA_LIBS})
  target_include_directories(${CASE_NAME} PRIVATE ${DCA_INCLUDE_DIRS})
  target_compile_definitions(${CASE_NAME} PRIVATE
    DCA_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")
  target_compile_definitions(${CASE_NAME} PRIVATE
                           ADIOS2_PERF_TEST_RANKS=${num_ranks})
endfunction()

adios_scaling_test(adios2 6)
adios_scaling_test(adios2 12)
adios_scaling_test(adios2 24)

endif()
