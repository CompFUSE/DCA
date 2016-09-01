# Builds main_dca.

if (DCA_BUILD_DCA)
  add_executable(main_dca main_dca.cpp)
  target_link_libraries(main_dca ${DCA_LIBRARIES})
  
  if (DCA_WITH_CUDA)
    target_link_libraries(main_dca
      DCA_GPU_routines
      ${DCA_GPU_LIBRARIES})
    cuda_add_cublas_to_target(main_dca)
  endif()

  target_include_directories(main_dca PRIVATE ${DCA_INCLUDES})
endif()
