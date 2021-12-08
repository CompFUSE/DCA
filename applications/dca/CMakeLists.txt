# Builds main_dca.

if (DCA_BUILD_DCA)
  add_executable(main_dca main_dca.cpp)
  target_include_directories(main_dca PRIVATE ${DCA_INCLUDE_DIRS})
  target_link_libraries(main_dca PRIVATE ${DCA_LIBS})
  
  if (DCA_WITH_CUDA)
    dca_gpu_blas_link(main_dca)
    target_compile_definitions(main_dca PUBLIC DCA_HAVE_CUDA)
  endif()
endif()
