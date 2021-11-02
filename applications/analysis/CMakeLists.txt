# Builds main_analysis.

if (DCA_BUILD_ANALYSIS)
  add_executable(main_analysis main_analysis.cpp)
  target_include_directories(main_analysis PRIVATE ${DCA_INCLUDE_DIRS})

  if (DCA_WITH_CUDA)
    target_link_libraries(main_dca PRIVATE CUDA::cublas)
  endif()

  if (DCA_WITH_HIP)
    target_compile_definitions(main_dca PUBLIC DCA_HAVE_HIP)
    target_link_libraries(main_dca PRIVATE hip::host roc::hipblas)
  endif()

  target_link_libraries(main_analysis ${DCA_LIBS})
endif()
