# test/performance/phys/tp_accumulation

if (DCA_WITH_TESTS_PERFORMANCE)
  add_executable(tp_accumulation_performance_test tp_accumulation_performance_test.cpp)
  target_link_libraries(tp_accumulation_performance_test PRIVATE ${DCA_LIBS})
  target_include_directories(tp_accumulation_performance_test PRIVATE  ${DCA_INCLUDE_DIRS})
  target_compile_definitions(tp_accumulation_performance_test PRIVATE
                             DCA_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")

  add_executable(sp_accumulation_performance_test sp_accumulation_performance_test.cpp)
  target_link_libraries(sp_accumulation_performance_test PRIVATE ${DCA_LIBS})
  target_include_directories(sp_accumulation_performance_test PRIVATE  ${DCA_INCLUDE_DIRS})
  target_compile_definitions(sp_accumulation_performance_test PRIVATE
    DCA_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")

  if (DCA_HAVE_GPU)
    dca_gpu_runtime_link(tp_accumulation_performance_test)
    dca_gpu_runtime_link(sp_accumulation_performance_test)
    if (DCA_HAVE_CUDA)
      target_link_libraries(tp_accumulation_performance_test PRIVATE nvToolsExt)
      target_link_libraries(sp_accumulation_performance_test PRIVATE nvToolsExt)
    endif()
  endif()
endif()
