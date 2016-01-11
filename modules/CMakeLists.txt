# modules/CMakeLists.txt

add_custom_target(
    modules_helper
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/cmake/modules_rebuild.cmake
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_library(modules ${CMAKE_CURRENT_BINARY_DIR}/modules.cpp)
target_include_directories(modules PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_dependencies(modules modules_helper)
