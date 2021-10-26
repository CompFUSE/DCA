
add_executable(ctint_walker_performance_test ctint_walker_performance_test.cpp)

target_link_libraries(ctint_walker_performance_test PRIVATE ${DCA_LIBS})
target_include_directories(ctint_walker_performance_test  PRIVATE
    ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
target_compile_definitions(ctint_walker_performance_test PRIVATE DCA_SOURCE_DIR="${PROJECT_SOURCE_DIR}")


if (DCA_HAVE_GPU)
  # target_compile_definitions(ctint_walker_performance_test PRIVATE DCA_HAVE_CUDA)
  dca_gpu_blas_link(ctint_walker_performance_test)
  target_link_libraries(ctint_walker_performance_test PRIVATE magma::magma mc_kernels)
endif ()
