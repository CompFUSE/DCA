# test/performance/phys/accumulation

if (DCA_WITH_TESTS_PERFORMANCE)
  add_executable(sp_accumulation_performance_test sp_accumulation_performance_test.cpp)
  target_link_libraries(sp_accumulation_performance_test PRIVATE ${DCA_LIBS})
  target_include_directories(sp_accumulation_performance_test PRIVATE  ${DCA_INCLUDE_DIRS})
  target_compile_definitions(sp_accumulation_performance_test PRIVATE
                             DCA_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")
endif()
