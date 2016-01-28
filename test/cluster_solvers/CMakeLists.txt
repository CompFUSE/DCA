################################################################################
# test/cluster_solvers/CMakeLists.txt
################################################################################

# CT-AUX pthreaded GPU
if (DCA_GPU_SUPPORT)
  add_executable(ctaux_pthreaded_gpu_test ${CMAKE_CURRENT_SOURCE_DIR}/ctaux_pthreaded_gpu_test.cpp)
  target_link_libraries(ctaux_pthreaded_gpu_test
    DCA_GPU_routines
    ${DCA_LIBRARIES}
    ${DCA_GPU_LIBRARIES}
    )
  cuda_add_cublas_to_target(ctaux_pthreaded_gpu_test)

else()
  add_executable(ctaux_pthreaded_gpu_test ${CMAKE_CURRENT_SOURCE_DIR}/ctaux_pthreaded_gpu_test.cpp)
  target_link_libraries(ctaux_pthreaded_gpu_test ${DCA_LIBRARIES})
endif()

target_include_directories(ctaux_pthreaded_gpu_test PRIVATE ${DCA_INCLUDES})

# CT-AUX pthreaded CPU
add_executable(ctaux_pthreaded_cpu_test ${CMAKE_CURRENT_SOURCE_DIR}/ctaux_pthreaded_cpu_test.cpp)
target_link_libraries(ctaux_pthreaded_cpu_test ${DCA_LIBRARIES})
target_include_directories(ctaux_pthreaded_cpu_test PRIVATE ${DCA_INCLUDES})
