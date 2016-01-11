################################################################################
# include/CMakeLists.txt
################################################################################

configure_file(include_files.h.in "${CMAKE_CURRENT_BINARY_DIR}/include_files.h" @ONLY)
file(COPY type_definitions.h DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY enumerations.hpp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

