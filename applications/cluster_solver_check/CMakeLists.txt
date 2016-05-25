################################################################################
# applications/cluster_solver_check/CMakeLists.txt
################################################################################

configure_file(include_files.inc.in "${CMAKE_CURRENT_BINARY_DIR}/include_files.inc" @ONLY)

if (DCA_BUILD_CLUSTER_SOLVER_CHECK)
  message(STATUS "Generating cluster_solver_check.cpp.")
  configure_file(cluster_solver_check.cpp.in "${CMAKE_CURRENT_BINARY_DIR}/cluster_solver_check.cpp" @ONLY)
  
  if (DCA_WITH_CUDA)
    add_executable(cluster_solver_check ${CMAKE_CURRENT_BINARY_DIR}/cluster_solver_check.cpp)
    target_link_libraries(cluster_solver_check
      DCA_GPU_routines
      ${DCA_LIBRARIES}
      ${DCA_GPU_LIBRARIES}
      )
    cuda_add_cublas_to_target(cluster_solver_check)
    
  else()
    add_executable(cluster_solver_check ${CMAKE_CURRENT_BINARY_DIR}/cluster_solver_check.cpp)
    target_link_libraries(cluster_solver_check ${DCA_LIBRARIES})
  endif()
  
  target_include_directories(cluster_solver_check
    PRIVATE "${DCA_APPLICATION_FOLDER}"
    PRIVATE ${DCA_INCLUDES}
    )
endif()

add_subdirectory(test/CT-AUX/bilayer_lattice)
add_subdirectory(test/CT-AUX/square_lattice)
