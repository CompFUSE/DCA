include_directories(
    ${DCA++_SOURCE_DIR}
    ${DCA++_SOURCE_DIR}/applications
    ${DCA++_SOURCE_DIR}/src
    ${DCA_INCLUDES}
)
    
set(DCA_FUNCTION_LIBRARY_TESTS
  function_library_test.cpp
)

add_executable(function_library_test_1
  ${DCA_FUNCTION_LIBRARY_TESTS}
)

target_link_libraries(function_library_test_1 
  ${DCA_LIBRARIES}
  ${HDF5_LIBRARY}
  gtest
) 

add_test(
  NAME dca_function_library_test
  COMMAND "$<TARGET_FILE:function_library_test_1>"
)
