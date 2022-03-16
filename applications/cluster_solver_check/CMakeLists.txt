# Builds cluster_solver_check.

if (DCA_BUILD_CLUSTER_SOLVER_CHECK)
  add_executable(cluster_solver_check cluster_solver_check.cpp)
  target_include_directories(cluster_solver_check PRIVATE ${DCA_INCLUDE_DIRS})
  target_link_libraries(cluster_solver_check  PUBLIC ${DCA_LIBS} dca_hdf5 statistical_testing)

  if (DCA_HAVE_GPU)
    dca_gpu_runtime_link(cluster_solver_check)
    dca_gpu_blas_link(cluster_solver_check)
  endif()
endif()
