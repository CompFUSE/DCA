add_executable(ctaux_walker_performance_test ctaux_walker_performance_test.cpp)

target_link_libraries(ctaux_walker_performance_test PRIVATE ${DCA_LIBS})
target_include_directories(ctaux_walker_performance_test  PRIVATE
    ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
target_compile_definitions(ctaux_walker_performance_test PRIVATE DCA_SOURCE_DIR="${PROJECT_SOURCE_DIR}")


if (DCA_HAVE_GPU)
  dca_gpu_runtime_link(ctaux_walker_performance_test)
  dca_gpu_blas_link(ctaux_walker_performance_test)
  if (DCA_HAVE_CUDA)
    target_link_libraries(ctaux_walker_performance_test PRIVATE nvToolsExt)
  endif()
endif ()
