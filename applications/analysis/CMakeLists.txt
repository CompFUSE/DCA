################################################################################
# applications/analysis/CMakeLists.txt
################################################################################

configure_file(include_files.h.in "${CMAKE_CURRENT_BINARY_DIR}/include_files.h" @ONLY)

if (DCA_BUILD_ANALYSIS)
  message(STATUS "Generating main_analysis.cpp.")
  configure_file(main_analysis.cpp.in "${CMAKE_CURRENT_BINARY_DIR}/main_analysis.cpp" @ONLY)
  
  if (DCA_GPU_SUPPORT)
    cuda_add_executable(main_analysis ${CMAKE_CURRENT_BINARY_DIR}/main_analysis.cpp)
    target_link_libraries(main_analysis
      DCA_GPU_routines
      ${DCA_LIBRARIES}
      ${DCA_GPU_LIBRARIES}
      )
    cuda_add_cublas_to_target(main_analysis)
    
  else()
    add_executable(main_analysis ${CMAKE_CURRENT_BINARY_DIR}/main_analysis.cpp)
    target_link_libraries(main_analysis ${DCA_LIBRARIES})
  endif()
  
  target_include_directories(main_analysis
    PRIVATE ${DCA_INCLUDES}
    PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/..")
endif()

add_subdirectory(test)
