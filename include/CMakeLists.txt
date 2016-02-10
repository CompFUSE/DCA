################################################################################
# include/CMakeLists.txt
################################################################################

# Copy type_definitions.h only during the initial configuration, otherwise
# changes to the 'build copy' might be overwritten when reconfiguring.
if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/type_definitions.h")
  file(COPY type_definitions.h DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

configure_file(include_files.h.in "${CMAKE_CURRENT_BINARY_DIR}/include_files.h" @ONLY)
file(COPY enumerations.hpp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

