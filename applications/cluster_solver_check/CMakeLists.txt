# Builds cluster_solver_check.

if (DCA_BUILD_CLUSTER_SOLVER_CHECK)
  add_executable(cluster_solver_check cluster_solver_check.cpp)

  if (DCA_WITH_CUDA)
    target_link_libraries(cluster_solver_check
      DCA_GPU_routines
      ${DCA_LIBRARIES}
      ${DCA_GPU_LIBRARIES})
    cuda_add_cublas_to_target(cluster_solver_check)

  else()
    target_link_libraries(cluster_solver_check ${DCA_LIBRARIES})
  endif()

  target_include_directories(cluster_solver_check
    PRIVATE ${DCA_LATTICE_MODELS_FOLDER}
    PRIVATE ${DCA_INCLUDES})
endif()
