################################################################################
# examples/check_solver/CMakeLists.txt
################################################################################

if (DCA_GPU_SUPPORT)
  add_executable(check_solver ${CMAKE_CURRENT_SOURCE_DIR}/check_solver.cpp)
  target_link_libraries(check_solver
    DCA_GPU_routines
    ${DCA_LIBRARIES}
    ${DCA_GPU_LIBRARIES}
    )
  cuda_add_cublas_to_target(check_solver)

else()
  add_executable(check_solver ${CMAKE_CURRENT_SOURCE_DIR}/check_solver.cpp)
  target_link_libraries(check_solver ${DCA_LIBRARIES})
endif()

target_include_directories(check_solver PRIVATE ${DCA_INCLUDES})
