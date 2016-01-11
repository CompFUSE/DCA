################################################################################
# examples/dca/CMakeLists.txt
################################################################################

# DCA+ main (DCA+ loop)
message(STATUS "Generating main_dca.cpp.")
configure_file(main_dca.cpp.in "${CMAKE_CURRENT_BINARY_DIR}/main_dca.cpp" @ONLY)

if (DCA_GPU_SUPPORT)
  add_executable(main_dca ${CMAKE_CURRENT_BINARY_DIR}/main_dca.cpp)
  target_link_libraries(main_dca
    DCA_GPU_routines
    ${DCA_LIBRARIES}
    ${DCA_GPU_LIBRARIES}
    )
  cuda_add_cublas_to_target(main_dca)

else()
  add_executable(main_dca ${CMAKE_CURRENT_BINARY_DIR}/main_dca.cpp)
  target_link_libraries(main_dca ${DCA_LIBRARIES})
endif()

target_include_directories(main_dca PRIVATE ${DCA_INCLUDES})
