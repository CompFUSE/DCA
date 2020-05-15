# Builds cluster_solver_check.

if (DCA_BUILD_CLUSTER_SOLVER_CHECK)
  add_executable(cluster_solver_check cluster_solver_check.cpp)
  target_include_directories(cluster_solver_check PRIVATE ${DCA_INCLUDE_DIRS})
  target_link_libraries(cluster_solver_check  PRIVATE ${DCA_LIBS} statistical_testing)

  if (DCA_WITH_CUDA)
    cuda_add_cublas_to_target(cluster_solver_check)
  endif()
endif()
