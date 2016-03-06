# gitVersion/CMakeLists.txt

# Note pass real SRC dir in to script

add_custom_target(
    gitVersion_helper
    COMMAND ${CMAKE_COMMAND} -DSCRIPT_SRC_DIR=${PROJECT_SOURCE_DIR} -P ${CMAKE_SOURCE_DIR}/cmake/gitVersion_rebuild.cmake
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_library(gitVersion ${CMAKE_CURRENT_BINARY_DIR}/gitVersion.cpp)
target_include_directories(gitVersion PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_dependencies(gitVersion gitVersion_helper)
