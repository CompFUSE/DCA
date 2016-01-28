################################################################################
# test/cluster_solvers/CMakeLists.txt
################################################################################

if (DCA_GPU_SUPPORT)
  add_executable(ctaux_pthreaded_test ${CMAKE_CURRENT_SOURCE_DIR}/ctaux_pthreaded_test.cpp)
  target_link_libraries(ctaux_pthreaded_test
    DCA_GPU_routines
    ${DCA_LIBRARIES}
    ${DCA_GPU_LIBRARIES}
    )
  cuda_add_cublas_to_target(ctaux_pthreaded_test)

else()
  add_executable(ctaux_pthreaded_test ${CMAKE_CURRENT_SOURCE_DIR}/ctaux_pthreaded_test.cpp)
  target_link_libraries(ctaux_pthreaded_test ${DCA_LIBRARIES})
endif()

target_include_directories(ctaux_pthreaded_test PRIVATE ${DCA_INCLUDES})
