# gitVersion/CMakeLists.txt

add_custom_target(
    gitVersion_helper
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/cmake/gitVersion_rebuild.cmake
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_library(gitVersion ${CMAKE_CURRENT_BINARY_DIR}/gitVersion.cpp)
target_include_directories(gitVersion PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_dependencies(gitVersion gitVersion_helper)
