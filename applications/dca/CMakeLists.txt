# Builds main_dca.

if (DCA_BUILD_DCA)
  add_executable(main_dca main_dca.cpp)
  target_include_directories(main_dca PRIVATE ${DCA_INCLUDE_DIRS})
  
  if (DCA_WITH_CUDA)
    target_link_libraries(main_dca PRIVATE CUDA::cublas)
  endif()

  if (DCA_WITH_HIP)
    target_compile_definitions(main_dca PUBLIC DCA_HAVE_HIP)
    target_link_libraries(main_dca PRIVATE hip::host roc::hipblas)
  endif()

  target_link_libraries(main_dca PRIVATE signals ${DCA_LIBS})
endif()
