
add_executable(ctint_walker_performance_test ctint_walker_performance_test.cpp)

target_link_libraries(ctint_walker_performance_test PRIVATE ${DCA_LIBS})
target_include_directories(ctint_walker_performance_test  PRIVATE
    ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
target_compile_definitions(ctint_walker_performance_test PRIVATE DCA_SOURCE_DIR="${PROJECT_SOURCE_DIR}")


if (DCA_HAVE_CUDA)
  target_compile_definitions(ctint_walker_performance_test PRIVATE DCA_HAVE_CUDA)
  target_link_libraries(ctint_walker_performance_test PRIVATE ${MAGMA_LIBRARY};${CUDA_cublas_LIBRARY})
endif ()
