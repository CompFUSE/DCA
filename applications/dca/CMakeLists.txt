# Builds main_dca.

if (DCA_BUILD_DCA)
  add_executable(main_dca main_dca.cpp)

  if (DCA_WITH_CUDA)
    target_link_libraries(main_dca
      DCA_GPU_routines
      ${DCA_LIBRARIES}
      ${DCA_GPU_LIBRARIES})
    cuda_add_cublas_to_target(main_dca)
  else()
    target_link_libraries(main_dca ${DCA_LIBRARIES})
  endif()

  target_include_directories(main_dca
    PRIVATE ${DCA_LATTICE_MODELS_FOLDER}
    PRIVATE ${DCA_INCLUDES})
endif()
