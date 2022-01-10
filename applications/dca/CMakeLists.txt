# Builds main_dca.

if (DCA_BUILD_DCA)
  add_executable(main_dca main_dca.cpp)
  target_include_directories(main_dca PRIVATE ${DCA_INCLUDE_DIRS})

  if (DCA_HAVE_GPU)
    target_link_libraries(main_dca PRIVATE ${DCA_GPU_LIBS})
  endif()

  target_link_libraries(main_dca PUBLIC signals ${DCA_LIBS})
endif()
