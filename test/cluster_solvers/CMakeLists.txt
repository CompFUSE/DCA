################################################################################
# test/cluster_solvers/CMakeLists.txt
################################################################################

if (DCA_PTHREADS)

  if (DCA_GPU_SUPPORT)
    # CT-AUX pthreaded GPU
    add_executable(ctaux_pthreaded_gpu_test ${CMAKE_CURRENT_SOURCE_DIR}/ctaux_pthreaded_gpu_test.cpp)
    target_link_libraries(ctaux_pthreaded_gpu_test
      DCA_GPU_routines
      ${DCA_LIBRARIES}
      ${DCA_GPU_LIBRARIES}
      )
    cuda_add_cublas_to_target(ctaux_pthreaded_gpu_test)
    target_include_directories(ctaux_pthreaded_gpu_test PRIVATE ${DCA_INCLUDES})
  endif()

  # CT-AUX pthreaded CPU
  add_executable(ctaux_pthreaded_cpu_test ${CMAKE_CURRENT_SOURCE_DIR}/ctaux_pthreaded_cpu_test.cpp)
  target_link_libraries(ctaux_pthreaded_cpu_test ${DCA_LIBRARIES})
  target_include_directories(ctaux_pthreaded_cpu_test PRIVATE ${DCA_INCLUDES})

endif()

# CT-AUX CPU
add_executable(ctaux_cpu_test ${CMAKE_CURRENT_SOURCE_DIR}/ctaux_cpu_test.cpp)
target_link_libraries(ctaux_cpu_test ${DCA_LIBRARIES})
target_include_directories(ctaux_cpu_test PRIVATE ${DCA_INCLUDES})
