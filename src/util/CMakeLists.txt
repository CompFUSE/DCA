# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)

# Git version
# Note: Pass real source and binary dirs to script.
add_custom_target(git_version_helper
  COMMAND ${CMAKE_COMMAND} -DSCRIPT_SRC_DIR=${PROJECT_SOURCE_DIR}
  -DSCRIPT_BIN_DIR=${PROJECT_BINARY_DIR} -P ${PROJECT_SOURCE_DIR}/cmake/git_version_rebuild.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

add_library(git_version ${CMAKE_CURRENT_BINARY_DIR}/git_version.cpp)
add_dependencies(git_version git_version_helper)

# Modules
add_custom_target(modules_helper
  COMMAND ${CMAKE_COMMAND} -DSCRIPT_SRC_DIR=${PROJECT_SOURCE_DIR}
  -DSCRIPT_BIN_DIR=${PROJECT_BINARY_DIR} -P ${PROJECT_SOURCE_DIR}/cmake/modules_rebuild.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

add_library(modules ${CMAKE_CURRENT_BINARY_DIR}/modules.cpp)
add_dependencies(modules modules_helper)
