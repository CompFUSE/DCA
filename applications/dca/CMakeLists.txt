# Builds main_dca.

if (DCA_BUILD_DCA)
  add_executable(main_dca main_dca.cpp)
  target_include_directories(main_dca PRIVATE ${DCA_INCLUDE_DIRS})
  target_link_libraries(main_dca PRIVATE ${DCA_LIBS} function)
  
  if (DCA_WITH_CUDA)
    target_link_libraries(main_dca PRIVATE CUDA::cublas)
  endif()
  
  target_link_libraries(main_dca PRIVATE signals ${DCA_LIBS})
endif()
